using HDF5, Dates, Statistics, DelimitedFiles

include("./generator_functions.jl")

##### all parameters on the spatial and temporal domains are set below

# Set longitude and latitude limits and spacing for the grid
lonmin, lonspacing, lonmax = 15, 1, 60
latmin, latspacing, latmax = 30, 1, 75

# Choose the time bounds of interest: DateTime(Year,Month,Day,Hour,Minute,Second) [24 Hour Format]
start_date = DateTime(2003, 7, 26, 0, 0, 0)
end_date = DateTime(2003, 8, 6, 0, 0, 0)
time_step = Hour(6) # time step in hours (must be a multiple of 6)

##### finish setting parameters

# Set up an array containing the temporal range for 𝕄
date_range = start_date:time_step:end_date
num_time_steps = length(date_range)

# Set up the grid struct and a dictionary for it
println("Setting up the grid...")
d, grid = make_dict_grid(lonmin, lonmax, lonspacing, latmin, latmax, latspacing)

# Read in the data and calculate the diffusion parameters
println("Reading data and calculating parameters...")
lons_data, lats_data, u_data, v_data, ϵ, a_init = read_data_and_get_parameters(grid, date_range)

# The value of a computed above is an initial estimate. After some experimentation, this value of a was chosen to roughly match the leading spatial/temporal eigenvalues of the inflated generator.
a = 0.0032

# Create a generator at each discrete time instance and append each one to Gvec as you go along
Gvec = []

println("Creating time-slice generators...")
@showprogress for τ ∈ 1:num_time_steps
    
    # Generate the zonal and meridional velocity component interpolants
    Iplt_zonal, Iplt_meridional = get_linear_interpolant(lons_data, lats_data, u_data[:,:,τ], v_data[:,:,τ])

    # Define a velocity vector function from the interpolants
    F(x) = [Iplt_zonal(x[1], x[2]), Iplt_meridional(x[1], x[2])]

    # Create the generator at this time step and attach it to Gvec
    G = make_generator(d, grid, F, ϵ)
    push!(Gvec, G)
  
end

# Assemble the full inflated generator from each individual generator
println("Making inflated generator...")
@time 𝐆 = make_inflated_generator(Gvec, time_step, a)

println("Computing inflated eigenvalues...")
@time Λ, V = eigs(𝐆, which=:LR, nev=20, maxiter=100000)

println("Plotting slices...")
# Plot the spectrum
@time plot_spectrum(grid, Λ, V)

# Plot slices of a spatial eigenvector of your choice
println("Plotting eigenvector time slices...")
vector_index_to_plot = 10
time_slice_spacing = 4 # Plot every time_slice_spacing-th slice after start_date (time gap of time_slice_spacing*time_step)
moviefilename = "Movie of 10th inflated generator eigenvector for the East Block.gif"

@time plot_slices(real.(V), vector_index_to_plot, time_slice_spacing, grid, date_range, :RdBu, moviefilename)
# We have to use real() when producing plots for V or else an error will be thrown when attempting to plot the eigenvectors, even if imag(V[:,vecnum]) = zeros(size(V,1)), as V is a matrix of complex type.

# Calculate SEBA vectors using a collection of eigenvectors; only use real-valued spatial eigenvectors in SEBA
println("Computing SEBA vectors...")
seba_inds = [1, 2, 4, 5, 9, 10]
Σ, ℛ = SEBA(real.(V[:, seba_inds])) # Again, we must take the real part of V or an error will be thrown when running SEBA().
println("The respective SEBA vector minima are ", minimum(Σ, dims=1))

# Plot slices of an individual SEBA vector
println("Plotting SEBA vector time slices...")
seba_index_to_plot = 4
moviefilename = "Movie of 4th SEBA vector for the East Block.gif"

@time plot_slices(Σ, seba_index_to_plot, time_slice_spacing, grid, date_range, :Reds, moviefilename)

# Save the results to HDF5 and JLD2 files 
# Data to save: Vectors of lon/lat ranges (or the full grid struct in JLD2), date range vector, eigenvalues and eigenvectors of the inflated generator and SEBA vectors
println("Saving variables...")
name_save_file = "InfGen_Results_EuroBlock_East"
@time save_results(grid, date_range, Λ, V, Σ, name_save_file)