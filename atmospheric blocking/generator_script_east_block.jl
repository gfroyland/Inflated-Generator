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
lons_data, lats_data, u_data, v_data, ϵ, a = read_data_and_get_parameters(grid, date_range)

# Set a to this value to better match the leading temporal and spatial eigenvalues of the inflated generator
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
@time Λ, V = eigs(𝐆, which=:LR, nev=10, maxiter=100000)

println("Plotting slices...")
# Plot the spectrum and obtain the list of real valued spatial eigenvectors for SEBA
spectrumpicname = "./atmospheric blocking/Inflated Generator Eigenvalue Spectrum for the East Block.png"
@time real_spat_inds = plot_spectrum_and_get_real_spatial_eigs(grid, Λ, V, spectrumpicname)

# Calculate SEBA vectors using a collection of eigenvectors; only use real-valued spatial eigenvectors in SEBA
println("Computing SEBA vectors...")
Σ, ℛ = SEBA(real.(V[:, real_spat_inds])) # We must insert the real part of V into SEBA() or an error will be thrown, even if imag(V[:,k]) = zeros(size(V,1)), as V is a matrix of complex type.
println("The respective SEBA vector minima are ", minimum(Σ, dims=1))

# Plot slices of the SEBA vector showing the East Block
println("Plotting SEBA vector time slices...")
index_to_plot = 3 # The third SEBA vector illustrates this block
time_slice_spacing = 4
picfilename = "./atmospheric blocking/The East Block illustrated through SEBA vector $index_to_plot.png"
moviefilename = "./atmospheric blocking/Movie of the East Block illustrated through SEBA vector $index_to_plot.mp4"
@time plot_slices(Σ, index_to_plot, time_slice_spacing, grid, date_range, :Reds, picfilename, moviefilename)

# Save the results to HDF5 and JLD2 files 
# Data to save: Vectors of lon/lat ranges (or the full grid struct in JLD2), date range vector, time slice spacing for plots, eigenvalues and eigenvectors of the inflated generator and SEBA vectors
println("Saving variables...")
filename = "./atmospheric blocking/InfGen_Results_EuroBlock_East"
@time save_results(grid, date_range, time_slice_spacing, Λ, V, Σ, filename)