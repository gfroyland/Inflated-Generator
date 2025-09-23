# This script performs the full inflated generator method for the 3D RBC velocity system with the Arnoldi method for eigensolving the generator and SEBA vector computation performed on a GPU.
# Start by instantiating the Julia environment and loading in all necessary packages and functions for this method.

using LinearAlgebra, Interpolations, HCubature, SparseArrays, Statistics, HDF5, JLD2, DelimitedFiles, ArnoldiMethod, Plots, PColorPlot, CUDA

# Load in the generator functions and the SEBA function
include("generator_functions.jl")
include("SEBA.jl")

# Name the path to the folder containing the RBC inflated generator repository
pathname = "./src/gpu_code/aNew-Sep20-3/"

# Name the path to the folder containing the RBC velocity data
path_to_data = "./RBC Velocity Data/"

# Define the orientation of the RBC velocity (set this to 1 for +v (positive velocity) or -1 for -v (negative velocity))
v_orientation = 1 # 1 for +v, -1 for -v

# Define the time range for this system
start_time = 2001
final_time = 2103
Δt = 3

time_range = start_time:Δt:final_time
num_time_steps = length(time_range)

# Create a grid and indexing for the spatial domain [xmin,xmax]x[ymin,ymax]x[zmin,zmax]
println("Setting up the grid...")

ℓ = 1/16 # Set the box side lengths (should be 1/(2^n) where n ∈ ℤ, n ≥ 1)

xmin, Δx, xmax = -4, ℓ, 4 
ymin, Δy, ymax = -4, ℓ, 4
zmin, Δz, zmax = 0, ℓ, 1
d, grid = make_dict_grid(xmin, xmax, Δx, ymin, ymax, Δy, zmin, zmax, Δz)

# Set the spatial diffusion parameter ϵ and an initial heuristic for the temporal diffusion strength a

println("Reading data and calculating parameters...")
paramtxtfile = pathname * "InfGen Parameters.txt"
x_data, y_data, z_data, ϵ, a_init = read_data_and_get_parameters(grid, time_range, path_to_data, paramtxtfile)

# Choose a value for the temporal diffusion parameter a (you can start with the initial heuristic and increase it later on)
a = 9.25

println("Making time slices of generator...")
# Create the time averaged generator
Gvec = []
time_spatgens_total = 0
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
    global time_spatgens_total += time_spatgens
    push!(Gvec, G)
end

open(paramtxtfile,"a") do file
    println(file,"The time taken to construct all $(num_time_steps) spatial generators with multithreading enabled is ",time_spatgens_total," seconds")
end

open(paramtxtfile,"a") do file
    println(file,"The selected value for a is ",a)
end

println("Making inflated generator...")
𝐆 = make_inflated_generator(Gvec, Δt, a)

println("Computing inflated generator eigenvalues...")
num_of_Λ = 1000
tol = √eps() # The default tol for the Arnoldi method is used now, change it to any level you wish

Λ, V = eigensolve_inflated_generator(𝐆, num_of_Λ, tol, paramtxtfile)

# Classify the eigenvalues computed and plot the spectrum

spectrumpicname = pathname * "RBC Inflated Generator Spectrum.png"
real_spat_inds, temp_inds, comp_inds = classify_eigs_and_plot_spectrum(grid, Λ, V, tol, paramtxtfile, spectrumpicname)

# Make a movie and Figure file of the leading non-trivial real-valued spatial eigenvector V[:, real_spat_inds[2]] along the x-y Midplane

index_to_plot = real_spat_inds[2]
col_scheme = :RdBu
picfilename = pathname * "RBC Leading Spatial Eigenvector xy Midplane.png"
moviefilename = pathname * "RBC Leading Spatial Eigenvector xy Midplane.mp4"
plot_slices(real.(V), index_to_plot, grid, time_range, col_scheme, picfilename, moviefilename)

# Compute SEBA using as many of the available real-valued spatial eigenvectors as you wish

num_of_SEBA = min(30,length(real_spat_inds)) # Set the number of real-valued spatial eigenvectors to take
time_seba = @elapsed Σ, ℛ = SEBA(real.(V[:, real_spat_inds[1:num_of_SEBA]]))
println("The respective SEBA vector minima are ", minimum(Σ, dims=1))

open(paramtxtfile,"a") do file
    println(file,"The time taken to compute $(num_of_SEBA) SEBA vectors is ",time_seba," seconds")
    println(file,"The respective minima of the $(num_of_SEBA) SEBA vectors are: ",minimum(Σ, dims=1))
end

# Calculate the maxima of the SEBA vectors generated

Σ_max = maximum(Σ, dims=2)

# Interpolate the SEBA vectors over a finer resolution (if desired)

ℓ_int = ℓ/2
x_full, y_full, z_full, Σ_int, Σ_int_max = interpolate_vecs(Σ, grid, time_range, ℓ_int)

# Make a movie and Figure file of the SEBA maxima using your interpolated data along the x-y Midplane

index_to_plot = 1
col_scheme = :Reds
picfilename = pathname * "RBC SEBA Maxima xy Midplane.png"
moviefilename = pathname * "RBC SEBA Maxima xy Midplane.mp4"
plot_interpolated_slices(Σ_int_max, index_to_plot, x_full, y_full, z_full, time_range, col_scheme, picfilename, moviefilename)

println("Saving variables...")
genfilename = pathname * "RBC_InfGen_Spatial_Generator_Data.jld2"
resultsfilename = pathname * "RBC_InfGen_Eigenbasis_Results.h5"
sebafilename = pathname * "RBC_InfGen_SEBA_Data.h5"
@time save_results(grid, time_range, Gvec, Δt, genfilename, Λ, V, real_spat_inds, temp_inds, comp_inds, resultsfilename, Σ, Σ_max, x_full, y_full, z_full, Σ_int, Σ_int_max, sebafilename)

println("3D Inflated Generator Computations Complete!")
