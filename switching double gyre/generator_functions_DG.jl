using Plots, LinearAlgebra, QuadGK, SparseArrays, Arpack, Statistics, HDF5, JLD2

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
    # Set up x_range and y_range (the x and y coordinates (respectively) of the centre points of each box in our grid)
    x_range = x_min+Δ_x/2:Δ_x:x_max-Δ_x/2
    y_range = y_min+Δ_y/2:Δ_y:y_max-Δ_y/2
    x_range = round.(x_range, digits=6)
    y_range = round.(y_range, digits=6)
    # Build an array of tuples with each one representing the Cartesian coordinates of each box centre
    i = 0
    temparray = []
    for x ∈ x_range
        for y ∈ y_range
            i += 1
            push!(temparray, ([x, y], i))
        end
    end
    # Create a dictionary for the centres, obtain its keys and set up the struct
    d = Dict(temparray)
    centres = Tuple.(collect(keys(d)))
    grid = Grid(centres, x_range, y_range, x_min, x_max, y_min, y_max, Δ_x, Δ_y)
    return d, grid
end

"`dg_velocity(t, x)` returns the velocity for the switching Double Gyre system at the point x and time t"
function dg_velocity(t, x)

    # Define the time-dependent switching double gyre vector field F(t,x) see [Atnip/Froyland/Koltai, 2024]
    r(t) = (1 / 2) * (1 + tanh(10 * (t - (1 / 2))))
    α(t) = (1 - 2 * r(t)) / (3 * (r(t) - 2) * (r(t) + 1))
    β(t) = (2 - 9 * α(t)) / 3
    return [(-π / 2) * sin(π * (α(t) * x[1]^2 + β(t) * x[1])) * cos(π * x[2] / 2), (2 * α(t) * x[1] + β(t)) * cos(π * (α(t) * x[1]^2 + β(t) * x[1])) * sin(π * x[2] / 2)]

end

"`get_parameters(grid, T_range)` computes values for the inflated generator diffusion parameters ϵ and a"
function get_parameters(grid, T_range)
    
    # Calculate the median of the speeds in the Double Gyre system
    F_median = median(norm(dg_velocity(t, x)) for t ∈ T_range for x ∈ grid.centres)
    println("The median of the speeds is... $F_median")

    # Calculate the spatial diffusion parameter ϵ
    ϵ = √(0.1 * F_median * (grid.Δ_x))
    println("The calculated ϵ value is... $ϵ")

    # Calculate the temporal diffusion strength a
    L_max_x = (grid.x_max - grid.x_min)
    L_max_y = (grid.y_max - grid.y_min)
    
    a = (T_range[end]-T_range[1]) * √(1.1 * F_median * (grid.Δ_x)) / (max(L_max_x,L_max_y))
    println("The initial a value is... $a")

    return ϵ, a

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
    tol = 1e-2 
    intorder = 1
    #start looping over c (centres of cells)
    for c ∈ centres
        rightc = round.(c + δx, digits=6)
        if rightc ∈ keys(d)  #check that the box on the right exists
            #compute the entry for G corresponding to flux through the right face
            #Because of the integration from 0 to 1 instead of 0 to the length of the face, instead of
            #dividing by volume, I multiply by (Δ_y/volume), which is 1/Δ_x.
            #similarly in the other faces below
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

    #place negative row sums on the diagonal of G so that the row sum of G is now zero.
    G = G - spdiagm(vec(sum(spdiagm(1 ./ volume) * G * spdiagm(volume), dims=2)))

    #adjust G by a similarity transformation to ensure that the matrix has row sum 0
    #this ensures leading right evec is constant and leading left evec is a *density* rather than a measure on each box.
    #see Lemma 4.7 in FJK'13.
    G = spdiagm(1 ./ volume) * G * spdiagm(volume)

    return G

end

"`make_inflated_generator(Gvec, Δt, a)` constructs the inflated generator from the time-indexed vector of generator matrices `Gvec` with discrete temporal spacing `Δt` and temporal diffusion parameter `a`"
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

"`plot_spectrum(grid, Λ, V, spectrumpicname)` plots the length(Λ) leading eigenvalues of the inflated generator, distinguishing spatial eigenvalues from temporal ones using the means of the variances of the inflated generator eigenvectors. This function also returns a vector containing the indices of the first two real-valued spatial eigenvectors (including the first trivial eigenvector) for use later when calling SEBA."
function plot_spectrum_and_get_real_spatial_eigs(grid, Λ, V, spectrumpicname)

    # Number of spatial grid points
    N = length(grid.x_range) * length(grid.y_range)
    # Number of time slices
    T = Int(size(V)[1] / N)
    # Number of computed eigenvectors 
    K = size(V, 2)

    # Calculate Means of Eigenvector Variances 
    averagespatialvariance = [mean([var(V[(t-1)*N+1:t*N, k]) for t = 1:T]) for k = 1:K]

    # Distinguish spatial eigenvalues from temporal ones
    spat_inds = findall(x->x>1e-10,averagespatialvariance)
    temp_inds = findall(x->x<1e-10,averagespatialvariance)

    # Trivial Λ_1 should be plotted as a spatial eigenvalue, but averagespatialvariance[1] ≈ 0, so correct this before plotting
    popfirst!(temp_inds) 
    append!(spat_inds,1)

    # Find the first two real-valued spatial eigenvectors (including trivial V_1) for SEBA. We only require these two eigenvectors for our SEBA calculations in the Double Gyre example.
    real_spat_inds = intersect(findall(x->x>1e-10,averagespatialvariance),findall(x->abs(x)<1e-12,imag(Λ)))

    # Include 1 in real_spat_inds as well, sort the array so that 1 is listed first, then only retain the first two indices for this example.
    append!(real_spat_inds,1)
    real_spat_inds = sort(real_spat_inds)
    real_spat_inds = real_spat_inds[1:2]

    # Plot the spectrum
    scatter(Λ[spat_inds], label="Spatial Λ_k", shape=:circle, mc=:blue, title="$(length(Λ)) eigenvalues with largest real part, a = $a", xlabel="Re(Λ_k)", ylabel="Im(Λ_k)")
    scatter!(Λ[temp_inds], label="Temporal Λ_k", shape=:xcross, mc=:red, msw=4)
    xlabel!("Re(Λ_k)")
    display(ylabel!("Im(Λ_k)"))

    savefig(spectrumpicname)

    return real_spat_inds

end

"`plot_slices(V, index_to_plot, time_slice_spacing, grid, T_range, col_scheme, moviefilename)` plots every `time_slice_spacing`-th time slice of the spacetime vector from the `index_to_plot` column in the matrix of spacetime vectors `V` (can be eigenvectors or SEBA vectors) on the grid `grid` over the time steps in T_range. A colour scheme (col_scheme) should be chosen by the user. The animation of the vector slices over time will be saved to a file named `moviefilename.gif`, and the image of vector slices will be saved to `picfilename.png`."
function plot_slices(V, index_to_plot, time_slice_spacing, grid, T_range, col_scheme, picfilename, moviefilename)

    # Define the numbers of spatial grid points and time slices
    spacelength = length(grid.x_range) * length(grid.y_range)
    T = length(T_range)

    # If we're plotting V (generator eigenvectors), it should be ℓ^2-normalized. This step is not necessary for SEBA vectors. 
    if col_scheme == :RdBu
        V = stack(normalize.(eachcol(V))) * √(size(V, 1))
    end

    # If we're plotting SEBA vectors, all vector entries below 0 should be replaced with 0.
    if col_scheme == :Reds
        V[V .< 0] .= 0
    end

    # create a T-vector of time-slices (copies of space)
    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    # find a common colour range if we're plotting eigenvectors, otherwise fix col_lims to (0, 1) if we're plotting SEBA vectors.
    if col_scheme == :Reds
        col_lims = (0, 1)
    else
        col_lims = (minimum((V[:, index_to_plot])),maximum((V[:, index_to_plot])))
    end

    # create an animation of frames of the eigenvector
    anim = @animate for t = 1:time_slice_spacing:T
        tm = T_range[t]
        contourf(grid.x_range, grid.y_range, reshape(sliceV[t][:, index_to_plot], length(grid.y_range), length(grid.x_range)), clims=col_lims, c=col_scheme, xlabel="x", ylabel="y", title="t = $tm", linewidth=0, levels=100)
    end
    display(gif(anim, moviefilename, fps=8))

    # plot individual time frames
    fig = []
    for t = 1:time_slice_spacing:T
        tm = T_range[t]
        push!(fig, contourf(grid.x_range, grid.y_range, reshape(sliceV[t][:, index_to_plot], length(grid.y_range), length(grid.x_range)), clims=col_lims, c=col_scheme, title="t = $tm", linewidth=0, levels=100, xlim=(0, 3), ylim=(0, 2), aspectratio=1, legend=:none))
    end
    display(plot(fig..., layout=(3, 4)))
    savefig(picfilename)

end

"`save_results(grid, T_range, time_slice_spacing, Λ, V, Σ, filename)` saves relevant data and results from the inflated generator calculations to HDF5 and JLD2 files for subsequent use and analysis. Data saved: Grid ranges in x and y (or the entire grid struct in JLD2), the temporal range, inflated generator eigenvalues and eigenvectors; and SEBA vectors obtained from the eigenvectors."
function save_results(grid, T_range, time_slice_spacing, Λ, V, Σ, filename)
    
    # Save data to a JLD2 file for use in Julia
    filename_JLD2 = filename * ".jld2"
    jldsave(filename_JLD2; grid, T_range, Λ, V, Σ)

    # Save data to an HDF5 file for use in MATLAB and other programming languages which may not be able to process JLD2 files
    filename_HDF5 = filename * ".h5"
    file_ID = h5open(filename_HDF5, "w")

    file_ID["x_range"] = grid.x_range
    file_ID["y_range"] = grid.y_range

    # The collect() function should be used to save T_range in order to avoid an error being thrown
    file_ID["T_range"] = collect(T_range) 
    file_ID["time_slice_spacing"] = time_slice_spacing

    # Complex valued data cannot be saved to an HDF5 file, so the real and imaginary parts of the eigenvalues and eigenvectors must be split and saved separately
    file_ID["Eigvals_Real"] = real.(Λ)
    file_ID["Eigvals_Imag"] = imag.(Λ)

    file_ID["Eigvecs_Real"] = real.(V)
    file_ID["Eigvecs_Imag"] = imag.(V)

    file_ID["SEBA"] = Σ

    close(file_ID)

end