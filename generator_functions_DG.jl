using Plots, LinearAlgebra, QuadGK, SparseArrays, Arpack, Statistics

#create a data structure for the grid
struct Grid
    centres
    x_range
    y_range
    x_min
    x_max
    y_min
    y_max
    Δ_x
    Δ_y
end

#include the Sparse EigenBasis Approximation function
include("SEBA.jl")

"`make_dict_grid` creates a dictionary to set up our indexing of the grid and fill the grid struct"
function make_dict_grid(x_min, x_max, Δ_x, y_min, y_max, Δ_y)
    x_range = x_min+Δ_x/2:Δ_x:x_max-Δ_x/2
    y_range = y_min+Δ_y/2:Δ_y:y_max-Δ_y/2
    x_range = round.(x_range, digits=6)
    y_range = round.(y_range, digits=6)
    i = 0
    temparray = []
    for x ∈ x_range
        for y ∈ y_range
            i += 1
            push!(temparray, ([x, y], i))
        end
    end
    d = Dict(temparray)
    centres = Tuple.(collect(keys(d)))
    grid = Grid(centres, x_range, y_range, x_min, x_max, y_min, y_max, Δ_x, Δ_y)
    return d, grid
end

"`make_generator(d, grid, F, ϵ)` creates a generator on the grid given by the pair `d`, `grid` for the vector field `F` with spatial diffusion parameter `ϵ`."
function make_generator(d, grid, F, ϵ)

    #create list of box centres
    centres = collect(keys(d))

    #compute volumes of cells
    volume = zeros(length(centres))
    for c ∈ centres
        volume[d[c]] = grid.Δ_x * grid.Δ_y
    end

    #create incremental vectors to access adjacent grid cells
    δx = [grid.Δ_x, 0]
    δy = [0, grid.Δ_y]

    #create list of box centres
    #create an empty array G (for Generator) to hold the flux values. 
    x_length = length(grid.x_range)
    y_length = length(grid.y_range)
    Gdim = x_length * y_length
    G = spzeros(Gdim, Gdim)

    #set up normal vectors to each of the 4 faces and a parametric function for each face
    rightnormal = δx / norm(δx)  #normal vector pointing to the right...similarly make another 3, one for each direction
    rightface(c, t) = c + (δx / 2 + δy / 2) - δy * t   #parameterise right face of box with centre t by a parameter t ranging from 0 to 1
    leftnormal = -δx / norm(δx)
    leftface(c, t) = c + (-δx / 2 + δy / 2) - δy * t
    uppernormal = δy / norm(δy)
    upperface(c, t) = c + (δx / 2 + δy / 2) - δx * t
    lowernormal = -δy / norm(δy)
    lowerface(c, t) = c + (δx / 2 - δy / 2) - δx * t

    #construct the generator matrix G
    tol = 1e-2  #hard coded for now, but may need to be adapted to integral magnitudes
    intorder = 1
    #start looping over c (centres of cells)
    for c ∈ centres
        rightc = round.(c + δx, digits=6)
        if rightc ∈ keys(d)  #check that the box on the right exists
            #compute the entry for G corresponding to flux through the right face
            #Because of the integration from 0 to 1 instead of 0 to the length of the face, instead of
            #dividing by volume, I multiply by (Δ_y/volume), which is 1/Δ_x.
            #similarly in the other faces below
            #G[d[c], d[rightc]] = (quadgk(t -> max(F(rightface(c, t)) ⋅ rightnormal, 0) + ϵ, 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_x
            G[d[c], d[rightc]] = (quadgk(t -> max(F(rightface(c, t)) ⋅ rightnormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_x + (ϵ^2 / (2 * (grid.Δ_x)^2))
        end
        leftc = round.(c - δx, digits=6)
        if leftc ∈ keys(d)
            G[d[c], d[leftc]] = (quadgk(t -> max(F(leftface(c, t)) ⋅ leftnormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_x + (ϵ^2 / (2 * (grid.Δ_x)^2))
        end
        upperc = round.(c + δy, digits=6)
        if upperc ∈ keys(d)
            G[d[c], d[upperc]] = (quadgk(t -> max(F(upperface(c, t)) ⋅ uppernormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_y + (ϵ^2 / (2 * (grid.Δ_y)^2))
        end
        lowerc = round.(c - δy, digits=6)
        if lowerc ∈ keys(d)
            G[d[c], d[lowerc]] = (quadgk(t -> max(F(lowerface(c, t)) ⋅ lowernormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_y + (ϵ^2 / (2 * (grid.Δ_y)^2))
        end
    end
    #G = G - spdiagm(vec(sum(G, dims=2)))

    #place negative row sums on the diagonal of G so that the row sum of G is now zero.
    G = G - spdiagm(vec(sum(spdiagm(1 ./ volume) * G * spdiagm(volume), dims=2)))

    #adjust G by a similarity transformation to ensure that the matrix has row sum 0
    #this ensures leading right evec is constant and leading left evec is a *density* rather than a measure on each box.
    #see Lemma 4.7 in FJK'13.
    G = spdiagm(1 ./ volume) * G * spdiagm(volume)

    return G

end

"`make_inflated_generator(Gvec, Δt, a)` constructs theinflated generator from the time-indexed vector of generator matrices `Gvec` with discrete temporal spacing `Δt` and temporal diffusion parameter `a`"
function make_inflated_generator(Gvec, Δt, a)

    #create Gspat
    𝐆spat = blockdiag(Gvec...)

    #create Gtemp (which will be a^2 * 𝐋 below)
    T = length(Gvec)
    #create 1D Laplace finite-difference matrix
    L = Tridiagonal(ones(T - 1), -2 * ones(T), ones(T - 1))
    #adjust endpoint values to conserve mass (row/col sums = 0)
    L[1, 1] = -1
    L[T, T] = -1
    #create spacetime Laplace (aka Ltemp) with kronecker product
    𝐋 = kron(L/(Δt)^2, one(Gvec[1]))     #second argument is identity matrix of size Gvec[1]
    #additively combine Gspat and Gtemp to create inflated generator
    𝐆 = 𝐆spat + a^2 * 𝐋

    return 𝐆
end

#FROM HERE WE SHOULD DECIDE WHAT TO INCLUDE IN THE PUBLIC VERSION

"`plot_spectrum_meanvar(grid, Λ, V)` plots the length(Λ) leading eigenvalues of the spectrum of the inflated generator, along with the meanvariance plot to test each eigenvalue to see if it is spatial or temporal."
function plot_spectrum_meanvar(grid, Λ, V)

    #plot inflated spectrum
    display(scatter(Λ, title="$(length(Λ)) eigenvalues with largest real part, a = $a"))

    N = length(grid.x_range) * length(grid.y_range)
    T = Int(size(V)[1] / N)

    K = size(V, 2)

    meanvariance = [mean([var(V[(t-1)*N+1:t*N, k]) for t = 1:T]) for k = 1:K]

    # Visualise the results

    x_disp = 1:K
    display(scatter(x_disp, meanvariance)) # Remove all similar cases of Plots.etc

    return meanvariance
    # Plot spectrum with circles and crosses to distinguish spatial vs temporal eigvals in automated fashion
    #should be close to zero for temporal eigenfunctions (and the trivial eigenfunction) and far from zero for spatial eigenfunctions

end

"`plot_slices(V, grid, vecnum)` plots the spacetime eigenvector from the `vecnum` column in the matrix of spacetime eigenvectors `V` on the grid `grid`."
function plot_slices(V, grid, vecnum) # Add colourscheme as argument, add T_range as argument for time
# If max_seba, V is the maximum SEBA vector, let vecnum = 1
# vecnum supplied, should be real and spatial (ideally Λ_2)
    spacelength = length(grid.x_range) * length(grid.y_range)
    T = Int(size(V)[1] / spacelength)

    # ℓ^2-normalize each column of V and scale so that entries are close to ±1 
    V = stack(normalize.(eachcol(V))) * sqrt(size(V, 1))

    #create a T-vector of time-slices (copies of space)
    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    # find a common colour range
    cmax = maximum(abs.(real(V[:, vecnum])))

    # create an animation of frames of the eigenvector
    anim = @animate for t = 1:T
        tm = (t - 1) / (T - 1)
        contourf(grid.x_range, grid.y_range, reshape(real.(sliceV[t][:, vecnum]), length(grid.y_range), length(grid.x_range)), clims=(-cmax, cmax), c=:RdBu, xlabel="x", ylabel="y", title="t = $tm", linewidth=0, levels=100)
    end
    display(gif(anim, fps=10))

    # plot individual time frames
    fig = []
    for t = 1:T
        tm = (t - 1) / (T - 1)
        push!(fig, contourf(grid.x_range, grid.y_range, reshape(real.(sliceV[t][:, vecnum]), length(grid.y_range), length(grid.x_range)), clims=(-cmax, cmax), c=cgrad(:RdBu, rev=true), title="t = $tm", linewidth=0, levels=100, xlim=(0, 3), ylim=(0, 2), aspectratio=1, legend=:none))
    end
    display(plot(fig[1:2:end]..., layout=(3, 4)))

end

function save_results(grid, Λ, V, Σ, filename)
# Add comments to explain what this function does
# Save the time domain as well
# No manipulation of V
# Use both HDF5 and JLD2 to save (split real and imaginary parts in HDF5 only); jldsave for JLD2
# e.g. jldsave("name.jld2"; Λ, V, Σ, T)
# to load: load("name.jld2",Λ) dictionary object, variable names are keys.
    file_ID = h5open(filename, "w")

    file_ID["longitude"] = grid.x_range
    file_ID["latitude"] = grid.y_range

    Σ_max = maximum(Σ, dims=2)

    file_ID["IDL_Lambda_Real"] = real.(Λ)
    file_ID["IDL_Lambda_Imag"] = imag.(Λ)

    spacelength = length(grid.x_range) * length(grid.y_range)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    sliceΣ = [Σ[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    for τ ∈ 1:T

        for φ ∈ 1:size(V, 2)
            varname = "IDL_v_" * lpad(φ, 2, "0") * "_t" * lpad(τ, 2, "0")
            file_ID[varname] = real.(reshape(sliceV[τ][:, φ], length(grid.y_range), length(grid.x_range)))
        end
        for ψ ∈ 1:size(Σ, 2)
            varname = "IDL_SEBA_" * lpad(ψ, 2, "0") * "_t" * lpad(τ, 2, "0")
            file_ID[varname] = real.(reshape(sliceΣ[τ][:, ψ], length(grid.y_range), length(grid.x_range)))
        end

        varname = "IDL_SEBA_Max_t" * lpad(τ, 2, "0")
        file_ID[varname] = real.(reshape(Σ_max[(τ-1)*spacelength.+(1:spacelength)], length(grid.y_range), length(grid.x_range)))
    end

    close(file_ID)

end
