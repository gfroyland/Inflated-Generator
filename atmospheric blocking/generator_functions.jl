using Dates, Plots, LinearAlgebra, Interpolations, QuadGK, SparseArrays, Arpack, Statistics, ProgressMeter, HDF5, JLD2, DelimitedFiles

#create a data structure for the grid
struct Grid
    centres
    lonrange
    latrange
    lonmin
    lonmax
    latmin
    latmax
    lonspacing
    latspacing
end

include("./SEBA.jl")

#create a dictionary to do the indexing we want and a grid struct
function make_dict_grid(lonmin, lonmax, lonspacing, latmin, latmax, latspacing)
    # Define ranges of the centres of each box (in longitudinal and latitudinal coordinates)
    lonrange = lonmin+lonspacing/2:lonspacing:lonmax-lonspacing/2
    latrange = latmin+latspacing/2:latspacing:latmax-lonspacing/2
    lonrange = round.(lonrange, digits=6)
    latrange = round.(latrange, digits=6)
    # Set up an array of tuples, each one containing the coordinates of a box centre in the form (lon,lat)
    i = 0
    temparray = []
    for lon ∈ lonrange
        for lat ∈ latrange
            i += 1
            push!(temparray, ([lon, lat], i))
        end
    end
    # Set up the dictionary, obtain the indexing keys for the dictionary and load the struct for the grid
    d = Dict(temparray)
    centres = Tuple.(collect(keys(d)))
    grid = Grid(centres, lonrange, latrange, lonmin, lonmax, latmin, latmax, lonspacing, latspacing)
    return d, grid
end

# This function is used to read in the ECMWF velocity data used for our generator calculations over a spatial domain defined by grid and a time interval defined by date_range. While reading this data in, we also take the opportunity to calculate the diffusion parameters ϵ and a for the generator in this function.
function read_data_and_get_parameters(grid, date_range)

    # We use this to convert longitudinal/latitudinal measures of distance into units of metres, as the velocity data is given in units of metres per second and we will perform all generator calculations using units of metres for distance rather than degrees longitude/latitude.
    deg2metr = 40075000 / 360

    # Begin by reading in the longitudinal and latitudinal ranges of the velocity data
    name_of_file = "./atmospheric blocking/data/ERA5_atmos_6HR_" * string(year(date_range[1])) * lpad(month(date_range[1]), 2, "0") * lpad(day(date_range[1]), 2, "0") * "_" * lpad(hour(date_range[1]), 2, "0") * "00.h5"
    file_ID = h5open(name_of_file)

    # M_lons and M_lats are vectors containing the longitudinal 
    # and latitudinal coordinates (respectively) used to produce the grid over 
    # which the wind velocity data is defined.

    # M_lons is given in units of degrees longitude East (negative values indicate
    # degrees longitude West), while M_lats is given in units of degrees 
    # latitude North (South if negative).

    M_lons = read(file_ID, "/longitude")
    M_lats = read(file_ID, "/latitude")
    M_lats = reverse(M_lats)

    # The reason that reverse() is used above is because the latitudinal coordinate vector
    # supplied by ECMWF provides the latitudinal range in descending rather than ascending order.
    # This vector (and the u/v component matrices in the second dimension) therefore need 
    # to be reversed, or else the convenience constructors used when producing the linear
    # interpolants later on will (in the same fashion as MATLAB's griddedInterpolant) throw 
    # an error because the grid coordinates have not been supplied in ascending order.

    close(file_ID)
    
    # Locate the indices relevant to the spatial ranges for the domain of interest 𝕄 defined in grid.
    M_lons_index_range = findall(x->(x>=grid.lonmin)&&(x<=grid.lonmax),M_lons)
    M_lats_index_range = findall(x->(x>=grid.latmin)&&(x<=grid.latmax),M_lats)

    # Find the number of time steps over the full temporal range defined in 𝕄.
    num_time_steps = length(date_range)

    # u_data and v_data are the zonal and meridional components (respectively)
    # of the wind velocity vector data defined over the full grid constructed using
    # M_lons and M_lats. They are each given in units of metres per second. 

    u_data = zeros(length(M_lons),length(M_lats),num_time_steps)
    v_data = zeros(length(M_lons),length(M_lats),num_time_steps)

    for τ ∈ 1:num_time_steps

        # Locate the data file within this directory and open it to access the data
        name_of_file = "./atmospheric blocking/data/ERA5_atmos_6HR_" * string(year(date_range[τ])) * lpad(month(date_range[τ]), 2, "0") * lpad(day(date_range[τ]), 2, "0") * "_" * lpad(hour(date_range[τ]), 2, "0") * "00.h5"
        file_ID = h5open(name_of_file)
    
        # Read in and save wind velocity component data for the entire spatial domain, used to produce linear interpolants later.
        u_now = read(file_ID, "/ucomp")
        u_now = reverse(u_now, dims=2)
        u_data[:,:,τ] = u_now

        v_now = read(file_ID, "/vcomp")
        v_now = reverse(v_now, dims=2)
        v_data[:,:,τ] = v_now
        
        close(file_ID)
    
    end

    # After loading in all of the necessary velocity data, we calculate the diffusion parameters ϵ and a.
    # Calculate the median value of the box side lengths ℓ (in metres) over our spatial grid in 𝕄.
    centres = grid.centres
    num_centres = length(centres)

    ℓ = [zeros(num_centres) ; (grid.latspacing)*(deg2metr)*ones(num_centres)]

    for c ∈ 1:num_centres
        ℓ[c] = grid.lonspacing * cosd(centres[c][2]) * deg2metr
    end

    ℓ_median = median(ℓ) 
    println("The calculated ℓ_median is... $ℓ_median")

    # Calculate the median of the speeds over 𝕄.
    v̄ = median(sqrt.(u_data[M_lons_index_range,M_lats_index_range,:].^2 + v_data[M_lons_index_range,M_lats_index_range,:].^2))
    println("The median of the speeds is... $v̄")

    # Calculate the spatial diffusion parameter ϵ
    ϵ = √(0.1*v̄*ℓ_median)
    println("The calculated ϵ value is... $ϵ")

    # Calculate the temporal diffusion strength a
    L_max_lon = (grid.lonmax - grid.lonmin)*cosd(grid.latmin)*deg2metr
    L_max_lat = (grid.latmax - grid.latmin)*deg2metr

    a = ((date_range[end]-date_range[1])/Day(1))*√(1.1*v̄*ℓ_median)/(max(L_max_lon,L_max_lat)) 
    println("The initial heuristic for a is... $a")

    # Return the velocity data and longitude/latitude arrays needed to produce the interpolants, along with the diffusion parameters.
    return M_lons, M_lats, u_data, v_data, ϵ, a

end

# This function is used to produce linear interpolants for the zonal and meridional velocity components from the ECMWF data on a single time slice.
# We define these interpolants over the entire spatial range of the ECMWF data, not just the spatial range defined for 𝕄.
# A linear order of accuracy for our velocity interpolation was deemed to be sufficient after some experimentation.
# Time is not a variable for the below interpolants, we can easily define new interpolants at each time step as we build the generator later on.
function get_linear_interpolant(lons_data, lats_data, u_data, v_data) # Read in and use the data saved earlier in read_data_and_get_parameters

    # Here, lons_data represents longitude and lats_data represents latitude
    # The grid coordinates can be kept in units of degrees longitude/latitude for now, this will not affect the construction of the generator later.

    Iplt_z = linear_interpolation((lons_data, lats_data), u_data)
    Iplt_m = linear_interpolation((lons_data, lats_data), v_data)

    # Return interpolants for the zonal and meridional (respectively) velocity vector components on this particular time slice.
    return Iplt_z, Iplt_m
end

# This function is used to build a generator for the spatial domain defined by d and grid under a dynamic system defined by F and incorporating diffusion of level ϵ on a particular time slice.
function make_generator(d, grid, F, ϵ)

    # The strategy is to follow Froyland/Junge/Koltai (2013).  
    # Spatial coordinates are in degrees lon/lat (info contained in d and grid); distances will be converted to metres when producing the generator
    # Velocities (the output of F) are in m.s⁻¹
    # ϵ is a diffusion parameter, making a crude finite-difference diffusion on a 5-point stencil 
    
    #conversion factor from degrees to metres:  mean circumference of Earth at the equator is (roughly) 40075000m and there are 2π radians (360 degrees) in a circle, so to obtain metres from degrees at the equator we multiply by 40075000/360
    deg2metr = 40075000/360

    #create list of box centres
    centres = collect(keys(d))
    
    # Calculate the volume (or area as the spatial domain is two-dimensional) of each box in m^2.
    volume = zeros(length(centres))
    for c ∈ centres
        volume[d[c]] = grid.lonspacing * grid.latspacing * cosd(c[2]) * deg2metr^2
    end

    #create basic lon,lat increment vectors to access adjacent grid cells
    Δlon = [grid.lonspacing, 0]
    Δlat = [0, grid.latspacing]

    #create list of box centres
    #create an array G (for Generator) to hold the flux values. 
    #G is the main array we want to compute
    lonlength = length(grid.lonrange)
    latlength = length(grid.latrange)
    Gdim = lonlength * latlength
    G = spzeros(Gdim, Gdim)

    #set up normal vectors to each of the 4 faces and a parametric function for each face
    rightnormal = Δlon / norm(Δlon)  #normal vector pointing to the right...similarly make another 3, one for each direction
    rightface(c, t) = c + (Δlon / 2 + Δlat / 2) - Δlat * t   #parameterise right face of box with centre t by a parameter t ranging from 0 to 1
    leftnormal = -Δlon / norm(Δlon)
    leftface(c, t) = c + (-Δlon / 2 + Δlat / 2) - Δlat * t
    uppernormal = Δlat / norm(Δlat)
    upperface(c, t) = c + (Δlon / 2 + Δlat / 2) - Δlon * t
    lowernormal = -Δlat / norm(Δlat)
    lowerface(c, t) = c + (Δlon / 2 - Δlat / 2) - Δlon * t
    
    #construct the generator matrix G
    tol = 1e-2 
    intorder = 1
    #start looping over c (centres of cells)
    for c ∈ centres
        rightc = round.(c + Δlon, digits=6)
        if rightc ∈ keys(d)  #check that the box on the right exists
            #compute the entry for G corresponding to flux through the right face
            #for the additional diffusion term I use the standard 5-point stencil finite-difference approximation (where the 5th diagonal element is taken care of later by ensuring row sum is zero)
            G[d[c], d[rightc]] = (quadgk(t -> max(F(rightface(c, t)) ⋅ rightnormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / (grid.lonspacing * deg2metr * cosd(c[2])) + (ϵ^2 / (2 * (grid.lonspacing *deg2metr* cosd(c[2]))^2))
        end
        leftc = round.(c - Δlon, digits=6)
        if leftc ∈ keys(d)
            G[d[c], d[leftc]] = (quadgk(t -> max(F(leftface(c, t)) ⋅ leftnormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / (grid.lonspacing *deg2metr* cosd(c[2])) + (ϵ^2 / (2 * (grid.lonspacing *deg2metr* cosd(c[2]))^2))
        end
        upperc = round.(c + Δlat, digits=6)
        if upperc ∈ keys(d)
            G[d[c], d[upperc]] = (quadgk(t -> max(F(upperface(c, t)) ⋅ uppernormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / (grid.latspacing *deg2metr) + (ϵ^2 / (2 * (grid.latspacing*deg2metr)^2))
        end
        lowerc = round.(c - Δlat, digits=6)
        if lowerc ∈ keys(d)
            G[d[c], d[lowerc]] = (quadgk(t -> max(F(lowerface(c, t)) ⋅ lowernormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / (grid.latspacing *deg2metr) + (ϵ^2 / (2 * (grid.latspacing*deg2metr)^2))
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

# This function is used to make the inflated generator by applying temporal diffusion to the vector of generator matrices Gvec. The level of temporal diffusion is defined by a and a temporal spacing of time_step is used to produce the spacetime Laplace operator.
function make_inflated_generator(Gvec, time_step, a)

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
    𝐋 = kron(L/(time_step/Day(1))^2, one(Gvec[1]))    #second argument is identity matrix of size Gvec[1]
    #additively combine Gspat and Gtemp to create inflated generator
    𝐆 = 𝐆spat + a^2 * 𝐋

    return 𝐆
end

# This function is used to plot the length(Λ) leading eigenvalues of the inflated generator, distinguishing spatial eigenvalues from temporal ones using the means of the variances of the inflated generator eigenvectors V, and saves the plot to `spectrumpicname.png`. This function also returns a vector containing the indices of all real-valued spatial eigenvectors for use later when calling SEBA.
function plot_spectrum_and_get_real_spatial_eigs(grid, Λ, V, spectrumpicname)

    # Number of spatial grid points
    N = length(grid.lonrange) * length(grid.latrange)
    # Number of time slices
    T = Int(size(V)[1] / N)
    # Number of computed eigenvectors 
    K = size(V, 2)

    # Calculate Means of Eigenvector Variances 
    averagespatialvariance = [mean([var(V[(t-1)*N+1:t*N, k]) for t = 1:T]) for k = 1:K]

    # Plot the Spectrum, distinguishing spatial eigenvalues from temporal ones
    spat_inds = findall(x->x>1e-10,averagespatialvariance)
    real_spat_inds = intersect(findall(x->x>1e-10,averagespatialvariance),findall(x->abs(x)<1e-12,imag(Λ)))
    temp_inds = findall(x->x<1e-10,averagespatialvariance)
    
    # Trivial Λ_1 should be plotted as a spatial eigenvalue, but averagespatialvariance[1] ≈ 0, so correct this before plotting
    popfirst!(temp_inds) 
    append!(spat_inds,1)

    # Include 1 in real_spat_inds as well, and sort the array so that 1 is listed first
    append!(real_spat_inds,1)
    real_spat_inds = sort(real_spat_inds)

    # Plot the spectrum and return real_spat_inds for use in SEBA later
    scatter(Λ[spat_inds], label="Spatial Λ_k", shape=:circle, mc=:blue, title="$(length(Λ)) eigenvalues with largest real part, a = $a", xlabel="Re(Λ_k)", ylabel="Im(Λ_k)", size=(700,500))
    scatter!(Λ[temp_inds], label="Temporal Λ_k", shape=:xcross, mc=:red, msw=4)
    xlabel!("Re(Λ_k)")
    display(ylabel!("Im(Λ_k)"))

    savefig(spectrumpicname)

    return real_spat_inds

end

# This function plots every `time_slice_spacing`-th time slice of the spacetime vector from the `index_to_plot` column in the matrix of spacetime vectors `V` (can be eigenvectors or SEBA vectors) on the grid `grid` over the time steps in T_range. A colour scheme (col_scheme) should be chosen by the user. The animation of the vector slices over time will be saved to a file named `moviefilename.gif`, and the image of slices will be saved to `picfilename.png`.
function plot_slices(V, index_to_plot, time_slice_spacing, grid, date_range, col_scheme, titleforplots, picfilename, moviefilename)

    # Define the numbers of spatial grid points and time slices
    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = length(date_range)

    # If we're plotting SEBA vectors, all vector entries below 0 should be replaced with 0.
    if col_scheme == :Reds
        V[V .< 0] .= 0
    end

    #create a T-vector of time-slices (copies of space)
    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    # find a common colour range if we're plotting eigenvectors, otherwise fix col_lims to (0, 1) if we're plotting SEBA vectors.
    if col_scheme == :Reds
        col_lims = (0, 1)
    else
        col_lims = (minimum((V[:, index_to_plot])),maximum((V[:, index_to_plot])))
    end

    # create an animation of frames of the eigenvector
    anim = @animate for t = 1:time_slice_spacing:T
        title_now = titleforplots * ", " * Dates.format(date_range[t], "dd/mm/yy HH:MM")
        contourf(grid.lonrange, grid.latrange, reshape(sliceV[t][:, index_to_plot], length(grid.latrange), length(grid.lonrange)), clims=col_lims, c=col_scheme, xlabel="̊ E", ylabel="̊ N", title=title_now, linewidth=0, levels=100)
    end
    display(gif(anim, moviefilename, fps=8))

    # plot individual time frames
    fig = []
    for t = 1:time_slice_spacing:T
        title_now = Dates.format(date_range[t], "dd/mm HH:MM")
        push!(fig, contourf(grid.lonrange, grid.latrange, reshape(sliceV[t][:, index_to_plot], length(grid.latrange), length(grid.lonrange)), clims=col_lims, c=col_scheme, title=title_now, linewidth=0, levels=100, aspectratio=1, legend=:none))
    end
    display(plot(fig..., layout=(3, 4),size=(800,600),plot_title=titleforplots))
    savefig(picfilename)
    
end

# This function saves relevant data and results from the inflated generator calculations to HDF5 and JLD2 files for subsequent use and analysis. Data saved: Longitudinal and latitudinal grid ranges (or the entire grid struct in JLD2); the date range, inflated generator eigenvalues and eigenvectors; and SEBA vectors obtained from the eigenvectors.
function save_results(grid, date_range, time_slice_spacing, Λ, V, Σ, filename)
    
    # Save data to a JLD2 file for use in Julia
    filename_JLD2 = filename * ".jld2"
    jldsave(filename_JLD2; grid, date_range, time_slice_spacing, Λ, V, Σ)

    # Save data to an HDF5 file for use in MATLAB and other programming languages which may not be able to process JLD2 files
    filename_HDF5 = filename * ".h5"
    file_ID = h5open(filename_HDF5, "w")

    file_ID["lonrange"] = grid.lonrange
    file_ID["latrange"] = grid.latrange

    # Convert date_range into an array of date and time strings as DateTime objects cannot be saved to HDF5 files.
    date_range_str = String[]
    for i = 1:length(date_range)
        push!(date_range_str,Dates.format(date_range[i],"dd/mm/yy HH:MM"))
    end

    # The collect() function should be used to save date_range_str in order to avoid an error being thrown
    file_ID["date_range"] = collect(date_range_str) 
    file_ID["time_slice_spacing"] = time_slice_spacing

    # Complex valued data cannot be saved to an HDF5 file, so the real and imaginary parts of the eigenvalues and eigenvectors must be split and saved separately
    file_ID["Eigvals_Real"] = real.(Λ)
    file_ID["Eigvals_Imag"] = imag.(Λ)

    file_ID["Eigvecs_Real"] = real.(V)
    file_ID["Eigvecs_Imag"] = imag.(V)

    file_ID["SEBA"] = Σ

    close(file_ID)

end