using HDF5
using Dates
using Statistics
using DelimitedFiles

include("generator_functions_lininterp.jl")

# x and y arrays remain the same at each time step

# Choose the time bounds of interest using the next two lines
# DateTime(Year,Month,Day,Hour,Minute,Second) [24 Hour Format]
start_date = DateTime(2003, 7, 26, 0, 0, 0)
end_date = DateTime(2003, 8, 6, 0, 0, 0)
time_step = Hour(6) # time step in hours (must be a multiple of 6)

date_range = start_date:time_step:end_date
num_time_steps = length(date_range)

# conversion factor from degrees to metres:  mean circumference of Earth at the equator is 40075000m and there are 2π radians (360 degrees) in a circle, so to obtain metres from degrees at the equator we multiply by 40075000/360
deg2metr = 40075000/360 

# Before we begin, we calculate a suitable value of ϵ for the generator
# calculations below, using a median of the norm of the wind velocity vectors
# from our data set over the time interval of choice over all space and time.

# set up lon,lat grid objects now for a smaller domain
# lonmin, lonmax, latmin and latmax must all be contained within the original data domain bounds
# The above parameters are used to reduce the size of the longitude and latitude range vectors.
# δ_x/0.25 and δ_y/0.25 must both be whole numbers and factors of 4*(lonmax-lonmin) and 4*(latmax-latmin) respectively.

println("Setting up the grid...")

lonmin, lonspacing, lonmax = 15, 1, 60
latmin, latspacing, latmax = 30, 1, 75

d, grid = make_dict_grid(lonmin, lonmax, lonspacing, latmin, latmax, latspacing)

println("Calculating a suitable value for ϵ...")

centres = grid.centres

ℓ = [zeros(length(centres)) ; (grid.latspacing)*(deg2metr)*ones(length(centres))]

for c ∈ 1:length(centres)
    ℓ[c] = grid.lonspacing * cosd(centres[c][2]) * deg2metr
end

ℓ_median = median(ℓ) 
println("The calculated ℓ_median is... $ℓ_median")

# Calculate median of the speeds within 𝕄
# Read in longitude and latitude data from our velocity files, and use these to find appropriate index ranges pertaining to the spatial extent of 𝕄.

name_of_file = "ERA5_atmos_6Hourly_Summer2003/ERA5_atmos_6HR_" * string(year(start_date)) * lpad(month(start_date), 2, "0") * lpad(day(start_date), 2, "0") * "_" * lpad(hour(start_date), 2, "0") * "00.h5"
file_ID = h5open(name_of_file)

M_lons = read(file_ID, "/longitude")
M_lats = read(file_ID, "/latitude")
M_lats = reverse(M_lats)

M_lons_index_range = findall(x->(x>=lonmin)&&(x<=lonmax),M_lons)
M_lats_index_range = findall(x->(x>=latmin)&&(x<=latmax),M_lats)

close(file_ID)

# This 3D array stores the observed velocity data over all available spatial and temporal data points within 𝕄
speeds_over_𝕄 = zeros(length(M_lons_index_range)*length(M_lats_index_range),length(date_range))

for τ ∈ 1:length(date_range)

    name_of_file = "ERA5_atmos_6Hourly_Summer2003/ERA5_atmos_6HR_" * string(year(date_range[τ])) * lpad(month(date_range[τ]), 2, "0") * lpad(day(date_range[τ]), 2, "0") * "_" * lpad(hour(date_range[τ]), 2, "0") * "00.h5"
    file_ID = h5open(name_of_file)

    u_data = read(file_ID, "/ucomp")
    u_data = reverse(u_data, dims=2)

    v_data = read(file_ID, "/vcomp")
    v_data = reverse(v_data, dims=2)

    # speeds_now contains all of the wind speeds calculated for this particular time instance...
    speeds_now = sqrt.(u_data[M_lons_index_range,M_lats_index_range].^2 .+ v_data[M_lons_index_range,M_lats_index_range].^2)
    
    # ... which are then fed into the full array used to compute the median.
    global speeds_over_𝕄[:,τ] = speeds_now[:]
    
    close(file_ID)

end

v̄ = median(speeds_over_𝕄)
println("The median of the speeds is... $v̄")
ϵ = sqrt(0.1*v̄*ℓ_median)
println("The calculated ϵ value is... $ϵ")

Gvec = []

@showprogress for datetime_now ∈ date_range

    println("Reading data...")

    name_of_file = "ERA5_atmos_6Hourly_Summer2003/ERA5_atmos_6HR_" * string(year(datetime_now)) * lpad(month(datetime_now), 2, "0") * lpad(day(datetime_now), 2, "0") * "_" * lpad(hour(datetime_now), 2, "0") * "00.h5"
    file_ID = h5open(name_of_file)

    @time begin

        # lons_data and lats_data are vectors containing the ranges of longitudinal 
        # and latitudinal coordinates (respectively) used to produce the grid over 
        # which the wind velocity data is defined.

        # lons_data is given in units of degrees longitude East (negative indicates
        # degrees longitude West), while lats_data is given in units of degrees 
        # latitude North (South if negative, we won't be dealing with this case yet).

        lons_data = read(file_ID, "/longitude")
        lats_data = read(file_ID, "/latitude")

        # u_data and v_data are the zonal and meridional components (respectively)
        # of the wind velocity vector data defined over the grid constructured using
        # lons_data and lats_data. They are each given in units of metres per second. 

        u_data = read(file_ID, "/ucomp")
        v_data = read(file_ID, "/vcomp")

    end

    println("Creating interpolant...")

    @time Iplt_zonal, Iplt_meridional = get_linear_interpolant(lons_data, lats_data, u_data, v_data)
    F(x) = [Iplt_zonal(x[1], x[2]), Iplt_meridional(x[1], x[2])]

    println("Creating generator...")

    @time G = make_generator(d, grid, F, ϵ)

    push!(Gvec, G)

    if datetime_now == start_date
        λ_now, v_now = eigs(G, which=:LR, nev=10, maxiter=100000)
        println(collect(λ_now))
    end

    close(file_ID)

end

L_max_lon = (grid.lonmax - grid.lonmin)*cosd(grid.latmin)*deg2metr
L_max_lat = (grid.latmax - grid.latmin)*deg2metr
a = ((end_date-start_date)/Day(1))*sqrt(1.1*v̄*ℓ_median)/(max(L_max_lon,L_max_lat)) 
println("The heuristic for a is... $a")

a = 0.0032

println("Making inflated generator...")
@time 𝐆 = make_inflated_generator(Gvec, time_step, a)

println("Computing inflated eigenvalues...")
@time Λ, V = eigs(𝐆, which=:LR, nev=20, maxiter=100000)
println(collect(Λ))

println("Plotting slices...")
@time plot_spectrum(grid, Λ, V)

vecnum = 10
time_hop = 4 # Plot every time_hopth slice after start_date (time gap of time_hop*time_step)
figlayout = (3, 4)
@time plot_slices(real.(V), vecnum, time_hop, grid, date_range, :RdBu, figlayout)

seba_inds = [1 ; 2 ; 4 ; 5 ; 9 ; 10]
Σ, ℛ = SEBA(real.(V[:, seba_inds]))
println("The respective SEBA vector minima are ", minimum(Σ, dims=1))

println("Plotting SEBA vector time slices...")

sebanum = 4
@time plot_slices(Σ, sebanum, time_hop, grid, date_range, :Reds, figlayout)

# Save the results to HDF5 and JLD2 files 
# Data to save: Vectors of lon/lat ranges, date range vector, eigenvalues and eigenvectors of the inflated generator and SEBA vectors
println("Saving variables...")
name_save_file = "InfGen_Results_EuroBlock_East"
@time save_results(grid, date_range, Λ, V, Σ, name_save_file)