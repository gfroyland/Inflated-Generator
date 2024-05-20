using HDF5, Dates, Statistics, DelimitedFiles

include("generator_functions.jl")

##### all parameters on the spatial and temporal domains are set below

# Set longitude and latitude limits and spacing for the grid
lonmin, lonspacing, lonmax = 15, 1, 60
latmin, latspacing, latmax = 30, 1, 75

# Choose the time bounds of interest: DateTime(Year,Month,Day,Hour,Minute,Second) [24 Hour Format]
start_date = DateTime(2003, 7, 26, 0, 0, 0)
end_date = DateTime(2003, 8, 6, 0, 0, 0)
time_step = Hour(6) # time step in hours (must be a multiple of 6)

##### finish setting parameters

date_range = start_date:time_step:end_date
num_time_steps = length(date_range)

println("Setting up the grid...")
d, grid = make_dict_grid(lonmin, lonmax, lonspacing, latmin, latmax, latspacing)

println("Reading data...")
u_data_over_𝕄, v_data_over_𝕄, ϵ, a = read_velocity_data(grid)

#L_max_lon = (grid.lonmax - grid.lonmin)*cosd(grid.latmin)*deg2metr
#L_max_lat = (grid.latmax - grid.latmin)*deg2metr
#a = ((end_date-start_date)/Day(1))*√(1.1*v̄*ℓ_median)/(max(L_max_lon,L_max_lat)) 
#println("The estimate for a is... $a")

Gvec = []
date_range = start_date:time_step:end_date

println("Creating time-slice generators...")
@showprogress for datetime_now ∈ date_range
    
    Iplt_zonal, Iplt_meridional = get_linear_interpolant(lons_data, lats_data, u_data, v_data)

    F(x) = [Iplt_zonal(x[1], x[2]), Iplt_meridional(x[1], x[2])]

    G = make_generator(d, grid, F, ϵ)
    push!(Gvec, G)
  
end

println("Making inflated generator...")
@time 𝐆 = make_inflated_generator(Gvec, time_step, a)

println("Computing inflated eigenvalues...")
@time Λ, V = eigs(𝐆, which=:LR, nev=20, maxiter=100000)
println(collect(Λ))

println("Plotting slices...")
@time plot_spectrum(grid, Λ, V)

vector_index_to_plot = 10
time_slice_spacing = 4 # Plot every time_hopth slice after start_date (time gap of time_hop*time_step)
figure_layout = (3, 4)
@time plot_slices(real.(V), vector_index_to_plot, time_slice_spacing, grid, date_range, :RdBu, figure_layout)

seba_inds = [1 ; 2 ; 4 ; 5 ; 9 ; 10]
Σ, ℛ = SEBA(real.(V[:, seba_inds]))
println("The respective SEBA vector minima are ", minimum(Σ, dims=1))

println("Plotting SEBA vector time slices...")
seba_index_to_plot = 4
@time plot_slices(Σ, seba_index_to_plot, time_slice_spacing, grid, date_range, :Reds, figure_layout)

# Save the results to HDF5 and JLD2 files 
# Data to save: Vectors of lon/lat ranges, date range vector, eigenvalues and eigenvectors of the inflated generator and SEBA vectors
println("Saving variables...")
name_save_file = "InfGen_Results_EuroBlock_East"
@time save_results(grid, date_range, Λ, V, Σ, name_save_file)