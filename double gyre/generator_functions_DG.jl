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
include("../SEBA.jl")

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

"`plot_spectrum(grid, Λ, V)` plots the length(Λ) leading eigenvalues of the spectrum of the inflated generator, distinguishing spatial eigenvalues from temporal ones using the meanvariance of eigenvectors test."
function plot_spectrum(grid, Λ, V)

    # Calculate Means of Eigenvector Variances 

    N = length(grid.x_range) * length(grid.y_range)
    T = Int(size(V)[1] / N)

    K = size(V, 2)

    meanvariance = [mean([var(V[(t-1)*N+1:t*N, k]) for t = 1:T]) for k = 1:K]

    # Plot the Spectrum, distinguishing spatial eigenvalues from temporal ones
    
    spat_inds = findall(x->x>1e-10,meanvariance)
    temp_inds = findall(x->x<1e-10,meanvariance)

    # Trivial Λ_1 should be plotted as a spatial eigenvalue, but meanvariance[1] ≈ 0, alleviorate this before plotting

    popfirst!(temp_inds) 
    append!(spat_inds,1)

    scatter(Λ[spat_inds], label="Spatial Λ_k", shape=:circle, mc=:blue, title="$(length(Λ)) eigenvalues with largest real part, a = $a", xlabel="Re(Λ_k)", ylabel="Im(Λ_k)")
    scatter!(Λ[temp_inds], label="Temporal Λ_k", shape=:xcross, mc=:red, msw=4)
    xlabel!("Re(Λ_k)")
    display(ylabel!("Im(Λ_k)"))

end

"`plot_slices(V, vecnum, grid, T_range, col_scheme)` plots the spacetime eigenvector from the `vecnum` column in the matrix of spacetime eigenvectors `V` on the grid `grid` over the time steps in T_range. A colour scheme (col_scheme) should be chosen by the user."
function plot_slices(V, vecnum, grid, T_range, col_scheme)

    spacelength = length(grid.x_range) * length(grid.y_range)
    T = length(T_range)

    # If we're plotting V, it should be ℓ^2-normalized before being passed in to this function. 
    # V = stack(normalize.(eachcol(V))) * sqrt(size(V, 1))

    #create a T-vector of time-slices (copies of space)
    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    # find a common colour range
    col_lims = (minimum((V[:, vecnum])),maximum((V[:, vecnum])))

    # create an animation of frames of the eigenvector
    anim = @animate for t = 1:T
        tm = T_range[t]
        contourf(grid.x_range, grid.y_range, reshape(sliceV[t][:, vecnum], length(grid.y_range), length(grid.x_range)), clims=col_lims, c=col_scheme, xlabel="x", ylabel="y", title="t = $tm", linewidth=0, levels=100)
    end
    display(gif(anim, fps=10))

    # plot individual time frames
    fig = []
    for t = 1:T
        tm = T_range[t]
        push!(fig, contourf(grid.x_range, grid.y_range, reshape(sliceV[t][:, vecnum], length(grid.y_range), length(grid.x_range)), clims=col_lims, c=col_scheme, title="t = $tm", linewidth=0, levels=100, xlim=(0, 3), ylim=(0, 2), aspectratio=1, legend=:none))
    end
    display(plot(fig[1:2:end]..., layout=(3, 4)))

end

"`save_results(grid, T_range, Λ, V, Σ, filename)` saves relevant data and results from the inflated generator calculations to HDF5 and JLD2 files for further use later. Data saved: Grid ranges in x and y, the temporal range, inflated generator eigenvalues and eigenvectors; and SEBA vectors obtained from the eigenvectors."
function save_results(grid, T_range, Λ, V, Σ, filename)
    
    # Save data to JLD2 file
    filename_JLD2 = filename * ".jld2"
    jldsave(filename_JLD2; grid, T_range, Λ, V, Σ)

    # Save data to HDF5 file
    filename_HDF5 = filename * ".h5"
    file_ID = h5open(filename_HDF5, "w")

    file_ID["x_range"] = grid.x_range
    file_ID["y_range"] = grid.y_range
    file_ID["T_range"] = collect(T_range) # The collect() function must be used or an error will be thrown

    file_ID["Eigvals_Real"] = real.(Λ)
    file_ID["Eigvals_Imag"] = imag.(Λ)

    file_ID["Eigvecs_Real"] = real.(V)
    file_ID["Eigvecs_Imag"] = imag.(V)

    file_ID["SEBA"] = Σ

    close(file_ID)

end
