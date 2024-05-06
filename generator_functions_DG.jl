using Plots, LinearAlgebra, QuadGK, SparseArrays, Arpack, Statistics, ProgressMeter

#create a data structure for the grid;  it might be expanded later to include e.g. G, spectrum, evecs, etc... or other items
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

include("SEBA.jl")

#create a dictionary to do the indexing we want and a grid struct
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

function make_generator(d, grid, F, ϵ)

    #create list of box centres
    centres = collect(keys(d))

    #compute volumes of cells
    volume = zeros(length(centres))
    for c ∈ centres
        volume[d[c]] = grid.Δ_x * grid.Δ_y
    end

    #create basic lon,lat increment vectors to access adjacent grid cells
    δx = [grid.Δ_x, 0]
    δy = [0, grid.Δ_y]

    #create list of box centres
    #create an array G (for Generator) to hold the flux values. 
    #G is the main array we want to compute
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
            #I add ϵ directly to the flux value, so it is also integrated from 0 to 1
            #Because of the integration from 0 to 1 instead of 0 to the length of the face, instead of
            #dividing by volume, I multiply by (Δ_y/volume), which is 1/Δ_x.
            #similarly in the other faces below
            #G[d[c], d[rightc]] = (quadgk(t -> max(F(rightface(c, t)) ⋅ rightnormal, 0) + ϵ, 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_x
            G[d[c], d[rightc]] = (quadgk(t -> max(F(rightface(c, t)) ⋅ rightnormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_x + (ϵ^2 / (2 * (grid.Δ_x)^2))
        end
        leftc = round.(c - δx, digits=6)
        if leftc ∈ keys(d)
            #G[d[c], d[leftc]] = (quadgk(t -> max(F(leftface(c, t)) ⋅ leftnormal, 0) + ϵ, 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_x
            G[d[c], d[leftc]] = (quadgk(t -> max(F(leftface(c, t)) ⋅ leftnormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_x + (ϵ^2 / (2 * (grid.Δ_x)^2))
        end
        upperc = round.(c + δy, digits=6)
        if upperc ∈ keys(d)
            #G[d[c], d[upperc]] = (quadgk(t -> max(F(upperface(c, t)) ⋅ uppernormal, 0) + ϵ, 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_y
            G[d[c], d[upperc]] = (quadgk(t -> max(F(upperface(c, t)) ⋅ uppernormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_y + (ϵ^2 / (2 * (grid.Δ_y)^2))
        end
        lowerc = round.(c - δy, digits=6)
        if lowerc ∈ keys(d)
            #G[d[c], d[lowerc]] = (quadgk(t -> max(F(lowerface(c, t)) ⋅ lowernormal, 0) + ϵ, 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_y
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

function make_dynamic_generator(Gvec)

    Gᴰ = mean(Gvec)
    return Gᴰ

end

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

function plot_9vecs_IDL(grid, Λ, V)

    spacelength = length(grid.x_range) * length(grid.y_range)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    cmax = maximum(abs.(real(V)), dims=1)

    for τ ∈ 1:T
    #anim = @animate for τ ∈ 1:T

        τ_m = (τ - 1) / (T - 1)
        P = []
        for k = 1:9
            push!(P, Plots.contourf(grid.x_range, grid.y_range, reshape(real.(sliceV[τ][:, k]), length(grid.y_range), length(grid.x_range)), clims=(-cmax[k], cmax[k]), c=:RdBu, xlabel="x", ylabel="y", title="v_$k, t = $τ_m", linewidth=0, levels=100))
        end
        display(Plots.plot(P..., layout=(3, 3)))
    end

    #gif(anim, "EuroBlock_0DayExt_9EigVecs_a3p5_30Jan24.gif", fps=2)
    #gif(anim, "DoubleGyre_9EigVecs_22Apr24_2.mp4", fps=2)

end

function plot_spatemp_IDL(grid, Λ, V)

    #plot inflated spectrum
    display(Plots.scatter(Λ, title="$(length(Λ)) eigenvalues with largest real part, a = $a"))

    N = length(grid.x_range) * length(grid.y_range)
    T = Int(size(V)[1] / N)

    K = size(V, 2)

    meanvariance = [mean([var(V[(t-1)*N+1:t*N, k]) for t = 1:T]) for k = 1:K]

    # Visualise the results

    x_disp = 1:K
    display(Plots.scatter(x_disp, meanvariance))

    return meanvariance

    #should be close to zero for temporal eigenfunctions (and the trivial eigenfunction) and far from zero for spatial eigenfunctions

end

function plot_SEBA_IDL(grid, Σ)

    # Prepare the slices of nine eigenvectors at each time step 
    spacelength = length(grid.x_range) * length(grid.y_range)
    T = Int(size(Σ)[1] / spacelength)

    sliceΣ = [Σ[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    # Plot All SEBA Evolutions

    for σ ∈ 1:size(Σ)[2]
        vidname = "DG_SEBA_3Vecs_" * lpad(σ, 1, "0") * "_a1p5_31Jan24.mp4"
        for t ∈ 1:T
        #anim = @animate for t ∈ 1:T

            τ_m = (t - 1) / (T - 1)
            display(Plots.contourf(grid.x_range, grid.y_range, reshape(sliceΣ[t][:, σ], length(grid.y_range), length(grid.x_range)), clims=(0, 1), c=:Reds, xlabel="x", ylabel="y", title="SEBA $σ, t = $τ_m", linewidth=0, levels=100))

        end
        #gif(anim, vidname, fps=2)
        println("SEBA Video $σ Complete.")
    end

    # Plot the Maxima of each SEBA Vector

    vidname = "DG_SEBAMax_3Vecs_a1p5_31Jan24.mp4"
    for 𝒯 ∈ 1:T
    #anim = @animate for 𝒯 ∈ 1:T

        τ_m = (𝒯 - 1) / (T - 1)
        display(Plots.contourf(grid.x_range, grid.y_range, reshape(maximum(sliceΣ[𝒯][:, :], dims=2), length(grid.y_range), length(grid.x_range)), clims=(0, 1), c=:Reds, xlabel="x", ylabel="y", title="Max SEBA, t = $τ_m", linewidth=0, levels=100))

    end
    #gif(anim, vidname, fps=2)
    println("Max SEBA Video Complete.")

end

function save_results(grid, Λ, V, Σ, filename)

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
