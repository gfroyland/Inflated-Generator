# This file contains the functions required to numerically implement the Inflated Generator method on the three-dimensional Rayleigh-Benard Convection (RBC) velocity data.
# This version does not involve utilising a GPU to boost the efficiency of the Arnoldi method when computing the eigenbasis of the inflated generator, for the benefit of those who do not have access to a computer or a HPC system with GPU cores.

using LinearAlgebra, Interpolations, HCubature, SparseArrays, Statistics, HDF5, JLD2, DelimitedFiles, ArnoldiMethod, Plots, PColorPlot

#create a data structure for the grid
struct Grid
    centres
    x_range
    y_range
    z_range
    x_min
    x_max
    y_min
    y_max
    z_min
    z_max
    Δx
    Δy
    Δz
end

#create a dictionary to do the indexing we want and a grid struct
function make_dict_grid(x_min, x_max, Δx, y_min, y_max, Δy, z_min, z_max, Δz)
    # Define ranges of the centres of each box
    x_range = x_min+Δx/2:Δx:x_max-Δx/2
    y_range = y_min+Δy/2:Δy:y_max-Δy/2
    z_range = z_min+Δz/2:Δz:z_max-Δz/2

    x_range = round.(x_range, digits=6)
    y_range = round.(y_range, digits=6)
    z_range = round.(z_range, digits=6)

    # Set up an array of triples (3-tuples), each one containing the coordinates of a box centre in the form (x, y, z)
    i = 0
    temparray = []
    for z ∈ z_range
        for x ∈ x_range
            for y ∈ y_range
                i += 1
                push!(temparray, ([x, y, z], i))
            end
        end
    end

    # Set up the dictionary, obtain the indexing keys for the dictionary and load the struct for the grid
    d = Dict(temparray)
    centres = Tuple.(collect(keys(d)))
    grid = Grid(centres, x_range, y_range, z_range, x_min, x_max, y_min, y_max, z_min, z_max, Δx, Δy, Δz)
    return d, grid
end

# This function is used to read in the spatial coordinate data for the RBC flow system and calculate the diffusion parameters ϵ and a for the inflated generator.
function read_data_and_get_parameters(grid, time_range, path_to_data, paramtxtfile)

    # Read in the x, y and z coordinates for the velocity data
    name_of_file = path_to_data * "snapData_" * string(time_range[1]) * ".h5"
    file_ID = h5open(name_of_file)

    # M_x, M_y and M_z are vectors containing the x, y and z 
    # coordinates (respectively) over the full RBC spatial domain M

    M_x = read(file_ID, "/X")
    M_y = read(file_ID, "/Y")
    M_z = read(file_ID, "/Z")

    close(file_ID)

    # Find the number of time steps over the full temporal range defined in our time-inflated domain 𝕄 = [time_range[1], time_range[end]] × M
    num_time_steps = length(time_range)
    
    # Calculate the diffusion parameters ϵ and a

    # To save time and memory, we have pre-calculated the median of the RBC velocity data over the entire spacetime domain and included it below:
    v̄ = 0.21951913908660203
    println("The median of the speeds is... $v̄")

    # Calculate the spatial diffusion parameter ϵ
    ϵ = √(0.1*v̄*(grid.Δx))
    println("The calculated ϵ value is... $ϵ")

    # Calculate the temporal diffusion strength a
    # Note: This is an initial heuristic for a, it is likely that you will need to fine tune this value after computing the eigenbasis of the inflated generator
    L_max_x = (grid.x_max - grid.x_min)
    L_max_y = (grid.y_max - grid.y_min)
    L_max_z = (grid.z_max - grid.z_min)

    a = sqrt(5/32)*(time_range[end]-time_range[1])*√(1.1*v̄*(grid.Δx))
    println("The initial heuristic for a is... $a")

    # Load the previously calculated parameters and other relevant information to the inflated generator calculations into your parameter text file
    open(paramtxtfile,"w") do file
        println(file,"A grid of ",length(grid.x_range)," × ",length(grid.y_range)," × ",length(grid.z_range)," boxes over the spatial domain [",grid.x_min,",",grid.x_max,"] × [",grid.y_min,",",grid.y_max,"] × [",grid.z_min,",",grid.z_max,"] has been defined at each of ",num_time_steps," time steps within the interval [",time_range[1],",",time_range[end],"] t_f")
        println(file,"The time step is ",time_range[2]-time_range[1]," t_f")
        println(file,"The side length of each box is ",grid.Δx)
        println(file,"The median of the speeds is ",v̄)
        println(file,"The calculated ϵ value is ",ϵ)
        println(file,"The initial heuristic for a is ",a)
    end

    # Return the spatial coordinate arrays for the velocity data needed to produce the interpolants, along with the diffusion parameters.
    # For the sake of saving time and conserving memory, we will load the velocity data in and produce velocity interpolants one time step at a time later.
    return M_x, M_y, M_z, ϵ, a

end

# This function is used to produce linear interpolants for the velocity components of the RBC data on a single time slice.
# We define these interpolants over the entire spatial range of the RBC data, not just the spatial range defined for 𝕄.
# A linear order of accuracy for our velocity interpolation was deemed to be sufficient after some experimentation.
# Time is not a variable for the below interpolants, we can easily define new interpolants at each time step as we build the generator later on.
# Read in the spatial coordinate data (x_data, y_data and z_data) and the velocity component data (u_data, v_data and w_data) for a particular time slice.
function get_linear_interpolant(x_data, y_data, z_data, u_data, v_data, w_data) 

    # We create one interpolant per velocity component (three overall)

    Iplt_u = linear_interpolation((x_data, y_data, z_data), u_data, extrapolation_bc=(Periodic(),Periodic(),Reflect()))
    Iplt_v = linear_interpolation((x_data, y_data, z_data), v_data, extrapolation_bc=(Periodic(),Periodic(),Reflect()))
    Iplt_w = linear_interpolation((x_data, y_data, z_data), w_data, extrapolation_bc=(Periodic(),Periodic(),Reflect()))

    # Return the interpolants for the velocity vector components at a particular time slice.
    return Iplt_u, Iplt_v, Iplt_w
end

# This function is used to build a generator for the spatial domain defined by d and grid under a flow system defined by F and incorporating diffusion of level ϵ on a particular time slice.
# This function incorporates multithreading for computation of the flux integrals used to build the spatial generator, owing to the large number of boxes defined within our spatial domain and therefore a large number of integrals to compute.
function make_generator(d, grid, F, ϵ)

    # The strategy is to follow Froyland/Junge/Koltai (2013).  
    # ϵ is a diffusion parameter, making a crude finite-difference diffusion on a 5-point stencil 
    
    #create list of box centres
    centres = collect(keys(d))
    
    # Calculate the volume of each box.
    volume = zeros(length(centres))
    for c ∈ centres
        volume[d[c]] = grid.Δx * grid.Δy * grid.Δz
    end

    #create basic x,y,z increment vectors to access adjacent grid cells
    δx = [grid.Δx, 0, 0]
    δy = [0, grid.Δy, 0]
    δz = [0, 0, grid.Δz]

    # Create length vectors for the x and y dimensions of the domain for periodic boundary conditions (not needed in z where we have Neumann boundary conditions)
    Lx = [(grid.x_max-grid.x_min), 0, 0]
    Ly = [0, (grid.y_max-grid.y_min), 0]

    # To reduce computation time owing to the large number of boxes in the domain (and therefore the large number of flux integrals to compute), we employ the following multithreaded computation technique.
    # Start by calculating Gdim, the size of the spatial generator in each dimension
    x_length = length(grid.x_range)
    y_length = length(grid.y_range)
    z_length = length(grid.z_range)
    Gdim = x_length * y_length * z_length

    # Initialise vectors to store the row and column indices of the flux values and another vector for the flux values themselves
    G_row_inds = zeros(Int64, 6*Gdim)
    G_col_inds = zeros(Int64, 6*Gdim)
    G_vals = zeros(6*Gdim)

    #set up normal vectors to each of the 6 faces of each box and a parametric function for each face
    # We take the x axis as the horizontal axis going into and out of the page, the y axis as the horizontal axis running along the page from left to right, and the z axis as the vertical axis.
    outernormal = δx / norm(δx) #normal vector pointing "out of the page"...similarly make another 5, one for each direction
    outerface(c, s, t) = c + [0, s, t] + (δx / 2 - δy / 2 - δz / 2) # representing the face in the direction pointing out of the page
    innernormal = -δx / norm(δx)
    innerface(c, s, t) = c + [0, s, t] + (-δx / 2 - δy / 2 - δz / 2)
    rightnormal = δy / norm(δy)  
    rightface(c, s, t) = c + [s, 0, t] + (-δx / 2 + δy / 2 - δz / 2)   
    leftnormal = -δy / norm(δy)
    leftface(c, s, t) = c + [s, 0, t] + (-δx / 2 - δy / 2 - δz / 2)
    uppernormal = δz / norm(δz)
    upperface(c, s, t) = c + [s, t, 0] + (-δx / 2 - δy / 2 + δz / 2)
    lowernormal = -δz / norm(δz)
    lowerface(c, s, t) = c + [s, t, 0] + (-δx / 2 - δy / 2 - δz / 2)

    #construct the generator matrix G properly with flux integrals
    tol = 1e-2
    #start looping over c (centres of cells) with multithreading enabled through Base.Threads
    Threads.@threads for ctr ∈ eachindex(centres)

        c = centres[ctr] # Obtain the coordinates for the centrepoint of the box we are currently working with
        inds = 6*(ctr-1) .+ collect(1:6) # Set the indices for the G_row_inds, G_col_inds and G_vals vectors (a maximum of six entries in each vector for each box)

        # Apply periodic BCs for x (outerc and innerc) and y (rightc and leftc), but not z (upperc and lowerc)

        outerc = [round.(c + δx, digits=6), round.(c + δx - Lx, digits=6)]
        if (~isempty(intersect(outerc,keys(d))))  #check that the box out of the page exists
            #compute the entry for G corresponding to flux through the face pointing out of the page
            #for the additional diffusion term I use the standard 5-point stencil finite-difference approximation (where the 5th diagonal element is taken care of later by ensuring row sum is zero)
            G_row_inds[inds[1]] = d[c]
            G_col_inds[inds[1]] = d[intersect(outerc,keys(d))[1]]
            G_vals[inds[1]] = (hcubature(t -> max(F(outerface(c, t[1], t[2])) ⋅ outernormal, 0), (0, 0), (grid.Δy, grid.Δz), rtol=tol, atol=tol)[1]) / (volume[d[intersect(outerc,keys(d))[1]]]) + (ϵ^2 / (2 * (grid.Δx)^2))
        end

        innerc = [round.(c - δx, digits=6), round.(c - δx + Lx, digits=6)]
        if (~isempty(intersect(innerc,keys(d))))
            G_row_inds[inds[2]] = d[c]
            G_col_inds[inds[2]] = d[intersect(innerc,keys(d))[1]]
            G_vals[inds[2]] = (hcubature(t -> max(F(innerface(c, t[1], t[2])) ⋅ innernormal, 0), (0, 0), (grid.Δy, grid.Δz), rtol=tol, atol=tol)[1]) / (volume[d[intersect(innerc,keys(d))[1]]]) + (ϵ^2 / (2 * (grid.Δx)^2))
        end

        rightc = [round.(c + δy, digits=6), round.(c + δy - Ly, digits=6)]
        if (~isempty(intersect(rightc,keys(d))))
            G_row_inds[inds[3]] = d[c]
            G_col_inds[inds[3]] = d[intersect(rightc,keys(d))[1]]
            G_vals[inds[3]] = (hcubature(t -> max(F(rightface(c, t[1], t[2])) ⋅ rightnormal, 0), (0, 0), (grid.Δx, grid.Δz), rtol=tol, atol=tol)[1]) / (volume[d[intersect(rightc,keys(d))[1]]]) + (ϵ^2 / (2 * (grid.Δy)^2))
        end

        leftc = [round.(c - δy, digits=6), round.(c - δy + Ly, digits=6)]
        if (~isempty(intersect(leftc,keys(d))))
            G_row_inds[inds[4]] = d[c]
            G_col_inds[inds[4]] = d[intersect(leftc,keys(d))[1]]
            G_vals[inds[4]] = (hcubature(t -> max(F(leftface(c, t[1], t[2])) ⋅ leftnormal, 0), (0, 0), (grid.Δx, grid.Δz), rtol=tol, atol=tol)[1]) / (volume[d[intersect(leftc,keys(d))[1]]]) + (ϵ^2 / (2 * (grid.Δy)^2))
        end

        upperc = round.(c + δz, digits=6)
        if upperc ∈ keys(d)
            G_row_inds[inds[5]] = d[c]
            G_col_inds[inds[5]] = d[upperc]
            G_vals[inds[5]] = (hcubature(t -> max(F(upperface(c, t[1], t[2])) ⋅ uppernormal, 0), (0, 0), (grid.Δx, grid.Δy), rtol=tol, atol=tol)[1]) / (volume[d[upperc]]) + (ϵ^2 / (2 * (grid.Δz)^2))
        end

        lowerc = round.(c - δz, digits=6)
        if lowerc ∈ keys(d)
            G_row_inds[inds[6]] = d[c]
            G_col_inds[inds[6]] = d[lowerc]
            G_vals[inds[6]] = (hcubature(t -> max(F(lowerface(c, t[1], t[2])) ⋅ lowernormal, 0), (0, 0), (grid.Δx, grid.Δy), rtol=tol, atol=tol)[1]) / (volume[d[lowerc]]) + (ϵ^2 / (2 * (grid.Δz)^2))
        end

    end

    # Remove all zero-valued elements from the index and value vectors
    G_row_inds = G_row_inds[G_row_inds .≠ 0]
    G_col_inds = G_col_inds[G_col_inds .≠ 0]
    G_vals = G_vals[G_vals .≠ 0]

    # Set up the generator matrix G using the index and value vectors (and make sure that G is a Gdim*Gdim matrix)
    G = sparse(G_row_inds, G_col_inds, G_vals, Gdim, Gdim)

    #place negative row sums on the diagonal of G so that the row sum of G is now zero.
    G = G - spdiagm(vec(sum(spdiagm(1 ./ volume) * G * spdiagm(volume), dims=2)))

    #adjust G by a similarity transformation to ensure that the matrix has row sum 0
    #this ensures leading right evec is constant and leading left evec is a *density* rather than a measure on each box.
    #see Lemma 4.7 in FJK'13.
    G = spdiagm(1 ./ volume) * G * spdiagm(volume)

    return G
end

# This function is used to build the spatial generator at every time step within time_range, by calling make_generator() num_time_steps times.
function make_generator_slices(d, grid, x_data, y_data, z_data, v_orientation, ϵ, time_range, path_to_data, paramtxtfile)

    Gvec = []
    time_spatgens_total = 0
    num_time_steps = length(time_range)

    for t ∈ 1:num_time_steps
        # Due to memory constraints, velocity data needs to be loaded in one time step at a time
        name_of_file = path_to_data * "snapData_" * string(time_range[t]) * ".h5"
        file_ID = h5open(name_of_file)

        u_now = read(file_ID, "/U")
        u_now = permutedims(u_now, [3, 2, 1])

        v_now = read(file_ID, "/V")
        v_now = permutedims(v_now, [3, 2, 1])

        w_now = read(file_ID, "/W")
        w_now = permutedims(w_now, [3, 2, 1])
        
        close(file_ID)

        # Generate the velocity component interpolants
        Iplt_u, Iplt_v, Iplt_w = get_linear_interpolant(x_data, y_data, z_data, u_now, v_now, w_now)

        # Define a velocity vector function from the interpolants
        F(x) = (v_orientation).*[Iplt_u(x[1], x[2], x[3]), Iplt_v(x[1], x[2], x[3]), Iplt_w(x[1], x[2], x[3])]

        # Create the generator at this time step and attach it to Gvec
        time_spatgens = @elapsed G = make_generator(d, grid, F, ϵ)
        time_spatgens_total += time_spatgens
        push!(Gvec, G)
    end

    open(paramtxtfile,"a") do file
        println(file,"The time taken to construct all $(num_time_steps) spatial generators with multithreading enabled is ",time_spatgens_total," seconds")
    end

    return Gvec

end

# This function is used to make the inflated generator by applying temporal diffusion to the vector of generator matrices Gvec. The level of temporal diffusion is defined by a and a temporal spacing of time_step is used to produce the spacetime Laplace operator.
function make_inflated_generator(Gvec, time_step, a)

    #create Gspat, a block diagonal matrix with the spatial generators at each time step along the diagonal
    𝐆spat = blockdiag(Gvec...)

    #create Gtemp (which will be a^2 * 𝐋 below)
    T = length(Gvec)
    #create 1D Laplace finite-difference matrix
    L = Tridiagonal(ones(T - 1), -2 * ones(T), ones(T - 1))
    #adjust endpoint values to conserve mass (row/col sums = 0)
    L[1, 1] = -1
    L[T, T] = -1
    #create spacetime Laplace (aka Ltemp) with kronecker product
    𝐋 = kron(L/(time_step)^2, one(Gvec[1]))    #second argument is identity matrix of size Gvec[1]
    #additively combine Gspat and Gtemp to create inflated generator
    𝐆 = 𝐆spat + a^2 * 𝐋

    return 𝐆
end

# This function uses the Arnoldi method to find the num_of_Λ eigenvalues and eigenvectors with largest real part of the inflated generator matrix 𝐆 for the RBC flow system.
# tol represents the tolerance for the eigenvalue/eigenvector residuals, and the computation time for the method is saved to paramtxtfile.
function eigensolve_inflated_generator(𝐆, num_of_Λ, tol, paramtxtfile)

    # Run the partial Schur decomposition component of the Arnoldi method
    time_decomp = @elapsed Decomp, History = partialschur(𝐆, nev=num_of_Λ, which=:LR, tol=tol)

    # Extract the eigenvalues and eigenvectors from Decomp
    time_eig = @elapsed Λ, V = partialeigen(Decomp)

    # The Λ vector and V matrix need to be reversed, so that the eigenvalues and eigenvectors are listed in the correct order
    Λ = reverse(Λ)
    V = reverse(V, dims=2)

    # Add a note on the computation time for the method to the parameter text file and return the eigenbasis
    open(paramtxtfile,"a") do file
        println(file,"The time taken to compute $num_of_Λ inflated generator eigenvalues using the Arnoldi Method on the CPU is ",(time_decomp+time_eig)," seconds")
    end

    return Λ, V

end

# This function is used to classify the length(Λ) leading eigenvalues of the inflated generator using the means of the variances of the inflated generator eigenvectors V. This function returns three vectors, with each one containing the indices of all real-valued spatial eigenvalues, temporal eigenvalues (if there are any) and complex valued eigenvalues respectively.
# We then plot the length(Λ) leading eigenvalues of the inflated generator, distinguishing spatial eigenvalues from temporal and complex ones using the index vectors generated, and return the index vectors for further use.
function classify_eigs_and_plot_spectrum(grid, Λ, V, tol, paramtxtfile, spectrumpicname)

    # Number of spatial grid points
    N = length(grid.x_range) * length(grid.y_range) * length(grid.z_range)
    # Number of time slices
    T = Int(size(V)[1] / N)
    # Number of computed eigenvectors 
    K = size(V, 2)

    # Calculate Means of Eigenvector Variances 
    averagespatialvariance = [mean([var(V[(t-1)*N+1:t*N, k]) for t = 1:T]) for k = 1:K]

    # Distinguish the spatial eigenvalues from the temporal and complex ones
    real_spat_inds = intersect(findall(x->x>tol,averagespatialvariance),findall(x->abs(x)<1e-12,imag(Λ)))
    temp_inds = findall(x->x<tol,averagespatialvariance)
    comp_inds = intersect(findall(x->x>tol,averagespatialvariance),findall(x->abs(x)>1e-12,imag(Λ)))

    # Trivial Λ_1 should be plotted as a spatial eigenvalue, but averagespatialvariance[1] ≈ 0, so check this before plotting
    # Add the lists of real-valued spatial and temporal eigenvalue indices to the parameter text file as you go along
    if (~isempty(temp_inds))
        
        if (abs(temp_inds[1]-1) < 1e-10)
            popfirst!(temp_inds) 
            append!(real_spat_inds,1)

            # Sort the real_spat_inds array so that 1 is listed first
            real_spat_inds = sort(real_spat_inds)
        end
        println("The temporal eigenvectors are: ",temp_inds)
        open(paramtxtfile,"a") do file
            println(file,"The temporal eigenvectors are: ",temp_inds)
        end

    else
        println("There are no temporal eigenvectors within the first $num_of_Λ for this value of a.")
        open(paramtxtfile,"a") do file
            println(file,"There are no temporal eigenvectors within the first $num_of_Λ for this value of a.")
        end
    end

    println("The usable (real-valued spatial) eigenvectors are: ",real_spat_inds)
    open(paramtxtfile,"a") do file
        println(file,"The usable (real-valued spatial) eigenvectors are: ",real_spat_inds)
    end

    # Plot the spectrum
    scatter(Λ[real_spat_inds], label="Spatial Λ_k", shape=:circle, mc=:blue, title="$(length(Λ)) eigenvalues with largest real part")
    if (~isempty(temp_inds))
        scatter!(Λ[temp_inds], label="Temporal Λ_k", shape=:xcross, mc=:red, msw=4)
    end
    scatter!(Λ[comp_inds], label="Complex Λ_k", shape=:square, mc=:black, msw=4, ma=0.2)
    scatter!(dpi=300)
    xlabel!("Re(Λ_k)")
    ylabel!("Im(Λ_k)")

    savefig(spectrumpicname)

    return real_spat_inds, temp_inds, comp_inds

end

# Plot eigenvector (or SEBA vector) no. index_to_plot along the x-y midplane
function plot_slices(V, index_to_plot, grid, time_range, col_scheme, picfilename, moviefilename)

    # Define the numbers of spatial grid points (in 2D and 3D) and time slices
    spacelength_2D = length(grid.x_range) * length(grid.y_range)
    spacelength_3D = length(grid.x_range) * length(grid.y_range) * length(grid.z_range)
    T = length(time_range)

    # If we're plotting SEBA vectors, all vector entries below 0 should be replaced with 0.
    if col_scheme == :Reds
        V[V .< 0] .= 0
    end

    # find a common colour range if we're plotting eigenvectors, otherwise fix col_lims to (0, 1) if we're plotting SEBA vectors.
    if col_scheme == :Reds
        col_lims = (0, 1)
    else
        col_lims = (minimum((V[:, index_to_plot])),maximum((V[:, index_to_plot])))
    end

    # Plot a movie and make a Figure for V (the vector in question needs to be re-sized first to make plotting easier)
    sliceV = [real.(V[(t-1)*spacelength_3D.+(1:spacelength_3D), :]) for t = 1:T]

    # Find the z index (or indices) corresponding to the x-y midplane

    if (mod(length(grid.z_range),2) == 0)
        ind_upper = trunc(Int64,ceil(length(grid.z_range)/2)+1)
        ind_lower = trunc(Int64,ceil(length(grid.z_range)/2))
    else
        ind_now = trunc(Int64,ceil(length(grid.z_range)/2))
    end

    # Make a movie of the eigenvector/SEBA vector of interest
    anim = @animate for t = 1:T
        title_now = "x-y Midplane, Time: $(time_range[t]) T_f"
        if (mod(length(grid.z_range),2) == 0)
            pcolor(grid.x_range, grid.y_range, ((reshape(sliceV[t][(ind_lower-1)*spacelength_2D.+(1:spacelength_2D), index_to_plot], length(grid.y_range), length(grid.x_range)) .+ reshape(sliceV[t][(ind_upper-1)*spacelength_2D.+(1:spacelength_2D), index_to_plot], length(grid.y_range), length(grid.x_range)))./2), interpolate=false, clims=col_lims, c=col_scheme, xlabel="x", ylabel="y", title=title_now, xlim=(-4, 4), ylim=(-4, 4), linewidth=0, levels=100, aspectratio=1, size=(600,500))
        else
            pcolor(grid.x_range, grid.y_range, reshape(sliceV[t][(ind_now-1)*spacelength_2D.+(1:spacelength_2D), index_to_plot], length(grid.y_range), length(grid.x_range)), interpolate=false, clims=col_lims, c=col_scheme, xlabel="x", ylabel="y", title=title_now, xlim=(-4, 4), ylim=(-4, 4), linewidth=0, levels=100, aspectratio=1, size=(600,500))
        end
    end
    gif(anim, moviefilename, fps=4)

    # Plot individual frames of the vector in a grid
    fig = []
    for t = 1:T
        title_now = "$(time_range[t]) T_f"
        if (mod(length(grid.z_range),2) == 0)
            push!(fig, pcolor(grid.x_range, grid.y_range, ((reshape(sliceV[t][(ind_lower-1)*spacelength_2D.+(1:spacelength_2D), index_to_plot], length(grid.y_range), length(grid.x_range)) .+ reshape(sliceV[t][(ind_upper-1)*spacelength_2D.+(1:spacelength_2D), index_to_plot], length(grid.y_range), length(grid.x_range)))./2), interpolate=false, clims=col_lims, c=col_scheme, title=title_now, xlim=(-4, 4), ylim=(-4, 4), linewidth=0, levels=100, aspectratio=1, legend=:none, size=(1000,800)))
        else
            push!(fig, pcolor(grid.x_range, grid.y_range, reshape(sliceV[t][(ind_now-1)*spacelength_2D.+(1:spacelength_2D), index_to_plot], length(grid.y_range), length(grid.x_range)), interpolate=false, clims=col_lims, c=col_scheme, title=title_now, xlim=(-4, 4), ylim=(-4, 4), linewidth=0, levels=100, aspectratio=1, legend=:none, size=(1000,800)))
        end
        
    end
    plot(fig..., layout=(5, 7), dpi=500)
    savefig(picfilename)

end

# Interpolate the vectors in V over a 3D spatial grid with box side length ℓ_int
function interpolate_vecs(V, grid, time_range, ℓ_int)

    # Define the new ranges of x, y and z coordinates over which we wish to interpolate
    x_full = (grid.x_min):ℓ_int:(grid.x_max)
    y_full = (grid.y_min):ℓ_int:(grid.y_max)
    z_full = (grid.z_min):ℓ_int:(grid.z_max)

    # Retrieve the numbers of time steps, grid points in (x,y) space and grid points in (x,y,z) space
    spacelength_2D = length(grid.x_range) * length(grid.y_range)
    spacelength_3D = length(grid.x_range) * length(grid.y_range) * length(grid.z_range)
    T = length(time_range)

    # Reshape the SEBA vectors to make interpolation a little easier
    sliceV = [V[(t-1)*spacelength_3D.+(1:spacelength_3D), :] for t = 1:T]

    # Perform linear interpolation of V with periodic BCs in x and y and reflective BCs in z
    V_int = zeros(length(x_full),length(y_full),length(z_full),T,size(V,2))

    for σ = 1:size(V,2)
        for t = 1:T
            v_now = zeros(length(grid.x_range),length(grid.y_range),length(grid.z_range))
            for z = 1:length(grid.z_range)
                ind_first = ((z-1)*spacelength_2D)+1;
                ind_last = (z*spacelength_2D);
                v_xy = sliceV[t][ind_first:ind_last, σ]
                v_now[:,:,z] = transpose(reshape(v_xy, length(grid.y_range), length(grid.x_range)))
            end
            v_fun = linear_interpolation((grid.x_range, grid.y_range, grid.z_range), v_now, extrapolation_bc = (Periodic(),Periodic(),Reflect()))
            V_int[:,:,:,t,σ] = v_fun(x_full,y_full,z_full)
        end
    end

    # Compute the maxima of the interpolated SEBA vectors
    V_int_max = maximum(V_int, dims=5)
    V_int_max = V_int_max[:,:,:,:,1] # Let Σ_int_max be a 4D array, not a 5D array with length 1 in the fifth dimension

    return x_full, y_full, z_full, V_int, V_int_max

end

# Plot eigenvector (or SEBA vector) no. index_to_plot along the x-y midplane
function plot_interpolated_slices(V, index_to_plot, x_full, y_full, z_full, time_range, col_scheme, picfilename, moviefilename)

    # Define the number of time slices
    T = length(time_range)

    # If we're plotting SEBA vectors, all vector entries below 0 should be replaced with 0.
    if col_scheme == :Reds
        V[V .< 0] .= 0
    end

    # find a common colour range if we're plotting eigenvectors, otherwise fix col_lims to (0, 1) if we're plotting SEBA vectors.
    if col_scheme == :Reds
        col_lims = (0, 1)
    else
        col_lims = (minimum((V[:,:,:,:,index_to_plot])),maximum((V[:,:,:,:,index_to_plot])))
    end

    # Find the z index (or indices) corresponding to the x-y midplane

    if (mod(length(z_full),2) == 0)
        ind_upper = trunc(Int64,ceil(length(z_full)/2)+1)
        ind_lower = trunc(Int64,ceil(length(z_full)/2))
    else
        ind_now = trunc(Int64,ceil(length(z_full)/2))
    end

    # Make a movie of the eigenvector/SEBA vector of interest
    anim = @animate for t = 1:T
        if (mod(length(z_full),2) == 0)
            V_now = (V[:,:,ind_lower,t,index_to_plot] .+ V[:,:,ind_upper,t,index_to_plot])./2
        else
            V_now = V[:,:,ind_now,t,index_to_plot]
        end

        title_now = "x-y Midplane, Time: $(time_range[t]) T_f"
        pcolor(x_full, y_full, transpose(V_now), interpolate=false, clims=col_lims, c=col_scheme, xlabel="x", ylabel="y", title=title_now, xlim=(-4, 4), ylim=(-4, 4), linewidth=0, levels=100, aspectratio=1, size=(600,500))
    end
    gif(anim, moviefilename, fps=4)

    # Plot individual frames of the vector in a grid
    fig = []
    for t = 1:T
        if (mod(length(z_full),2) == 0)
            V_now = (V[:,:,ind_lower,t,index_to_plot] .+ V[:,:,ind_upper,t,index_to_plot])./2
        else
            V_now = V[:,:,ind_now,t,index_to_plot]
        end

        title_now = "$(time_range[t]) T_f"
        push!(fig, pcolor(x_full, y_full, transpose(V_now), interpolate=false, clims=col_lims, c=col_scheme, title=title_now, xlim=(-4, 4), ylim=(-4, 4), linewidth=0, levels=100, aspectratio=1, legend=:none, size=(1000,800)))
    end
    plot(fig..., layout=(5, 7), dpi=500)
    savefig(picfilename)

end

# Plot eigenvector (or SEBA vector) no. index_to_plot along the x-y floor, ceiling and midplane; as well as the x-z and y-z vertical midplanes
function plot_interpolated_slices_fiveplanes(V, index_to_plot, x_full, y_full, z_full, time_ind, col_scheme, picfilename)

    # If we're plotting SEBA vectors, all vector entries below 0 should be replaced with 0.
    if col_scheme == :Reds
        V[V .< 0] .= 0
    end

    # find a common colour range if we're plotting eigenvectors, otherwise fix col_lims to (0, 1) if we're plotting SEBA vectors.
    if col_scheme == :Reds
        col_lims = (0, 1)
    else
        col_lims = (minimum((V[:,:,:,:,index_to_plot])),maximum((V[:,:,:,:,index_to_plot])))
    end

    # Plot individual frames of the vector along the five restricted 2D planes
    fig = []
    fig_layout = @layout [a b c; d ; e]
    
    # Plane #1: The x-y floor (Bottom)

    ind_now = 3 # Plot the max. SEBA not at the precise "floor" of the domain, but rather at a distance of twice the mesh size of the interpolated data from the floor
    V_now = V[:,:,ind_now,time_ind,index_to_plot]
    title_now = "xy-Bottom"
    push!(fig, pcolor(x_full, y_full, transpose(V_now), interpolate=false, clims=col_lims, c=col_scheme, title=title_now, xlim=(-4, 4), ylim=(-4, 4), linewidth=0, levels=100, aspectratio=1, legend=:none))
    
    # Plane #2: The x-y midplane
    # Find the z index (or indices) corresponding to the x-y midplane

    if (mod(length(z_full),2) == 0)
        ind_upper = trunc(Int64,ceil(length(z_full)/2)+1)
        ind_lower = trunc(Int64,ceil(length(z_full)/2))
        V_now = (V[:,:,ind_lower,time_ind,index_to_plot] .+ V[:,:,ind_upper,time_ind,index_to_plot])./2
    else
        ind_now = trunc(Int64,ceil(length(z_full)/2))
        V_now = V[:,:,ind_now,time_ind,index_to_plot]
    end

    title_now = "xy-Midplane"
    push!(fig, pcolor(x_full, y_full, transpose(V_now), interpolate=false, clims=col_lims, c=col_scheme, title=title_now, xlim=(-4, 4), ylim=(-4, 4), linewidth=0, levels=100, aspectratio=1, legend=:none))
    
    # Plane #3: The x-y ceiling (Top)

    ind_now = length(z_full)-2 # Plot the max. SEBA not at the precise "ceiling" of the domain, but rather at a distance of twice the mesh size of the interpolated data from the ceiling
    V_now = V[:,:,ind_now,time_ind,index_to_plot]
    title_now = "xy-Top"
    push!(fig, pcolor(x_full, y_full, transpose(V_now), interpolate=false, clims=col_lims, c=col_scheme, title=title_now, xlim=(-4, 4), ylim=(-4, 4), linewidth=0, levels=100, aspectratio=1, legend=:none))
    
    # Plane #4: The x-z midplane
    # Find the y index (or indices) corresponding to the x-z midplane

    V_perm = permutedims(V[:,:,:,time_ind,index_to_plot],[1 3 2])

    if (mod(length(y_full),2) == 0)
        ind_upper = trunc(Int64,ceil(length(y_full)/2)+1)
        ind_lower = trunc(Int64,ceil(length(y_full)/2))
        V_now = (V_perm[:,:,ind_lower] .+ V_perm[:,:,ind_upper])./2
    else
        ind_now = trunc(Int64,ceil(length(y_full)/2))
        V_now = V_perm[:,:,ind_now]
    end

    title_now = "xz-Midplane"
    push!(fig, pcolor(x_full, z_full, transpose(V_now), interpolate=false, clims=col_lims, c=col_scheme, title=title_now, xlim=(-4, 4), ylim=(0, 1), linewidth=0, levels=100, aspectratio=1, legend=:none))
    
    # Plane #5: The y-z midplane
    # Find the x index (or indices) corresponding to the y-z midplane

    V_perm = permutedims(V[:,:,:,time_ind,index_to_plot],[2 3 1])

    if (mod(length(x_full),2) == 0)
        ind_upper = trunc(Int64,ceil(length(x_full)/2)+1)
        ind_lower = trunc(Int64,ceil(length(x_full)/2))
        V_now = (V_perm[:,:,ind_lower] .+ V_perm[:,:,ind_upper])./2
    else
        ind_now = trunc(Int64,ceil(length(x_full)/2))
        V_now = V_perm[:,:,ind_now]
    end

    title_now = "yz-Midplane"
    push!(fig, pcolor(y_full, z_full, transpose(V_now), interpolate=false, clims=col_lims, c=col_scheme, title=title_now, xlim=(-4, 4), ylim=(0, 1), linewidth=0, levels=100, aspectratio=1, legend=:none))

    # Complete the full Figure using the above five plots
    plot(fig...; layout=fig_layout,size=(1000,1000),dpi=500)
    savefig(picfilename)

end

# This alternative version of save_results saves only the real-valued spatial eigenvectors to an HDF5 file, to avoid creating unworkably large results files (79 GB or more)
function save_results(grid, time_range, Gvec, Δt, genfilename, Λ, V, real_spat_inds, temp_inds, comp_inds, resultsfilename, Σ, Σ_max, x_full, y_full, z_full, Σ_int, Σ_int_max, sebafilename)
    
    # Three files to save: Spatial Generator (JLD2), Eigenbasis Results (HDF5), SEBA Results (HDF5)

    # File 1: Save the spatial generators to a JLD2 file so that the inflated generator calculations can be repeated for a new value of a without recalculating them
    jldsave(genfilename; grid, time_range, Gvec, Δt)

    # File 2: Save inflated generator eigenbasis data to an HDF5 file
    file_ID = h5open(resultsfilename, "w")

    # Save the spatial coordinate vectors
    file_ID["x_range"] = grid.x_range
    file_ID["y_range"] = grid.y_range
    file_ID["z_range"] = grid.z_range

    # The collect() function should be used to save time_range in order to avoid an error being thrown
    file_ID["time_range"] = collect(time_range)

    # Complex valued data cannot be saved to an HDF5 file, so the real and imaginary parts of the eigenvalues and eigenvectors must be split and saved separately
    file_ID["Eigvals_Real"] = real.(Λ)
    file_ID["Eigvals_Imag"] = imag.(Λ)

    # Only save the real-valued spatial eigenvectors, which for the 3D RBC flow likely won't even make up 10% of the overall eigenbasis sample
    file_ID["Eigvecs_Usable"] = real.(V[:, real_spat_inds])

    # Save the temp_inds, real_spat_inds and comp_inds vectors, so that a spectrum plot can still be reproduced
    file_ID["temp_inds"] = temp_inds
    file_ID["real_spat_inds"] = real_spat_inds
    file_ID["comp_inds"] = comp_inds

    close(file_ID)

    # File 3: Save SEBA data (original and interpolated) to an HDF5 file
    file_ID = h5open(sebafilename, "w")

    # Save the spatial coordinate vectors and the time range
    file_ID["x_range"] = grid.x_range
    file_ID["y_range"] = grid.y_range
    file_ID["z_range"] = grid.z_range
    file_ID["time_range"] = collect(time_range)

    # Save the original (non-interpolated) SEBA data
    file_ID["SEBA"] = Σ
    file_ID["SEBA_Max"] = Σ_max

    # Save the spatial coordinate vectors for the interpolated data
    file_ID["x_interp"] = collect(x_full)
    file_ID["y_interp"] = collect(y_full)
    file_ID["z_interp"] = collect(z_full)

    # Save the interpolated SEBA data
    file_ID["SEBA_Interp"] = Σ_int
    file_ID["SEBA_Interp_Max"] = Σ_int_max

    close(file_ID)

end