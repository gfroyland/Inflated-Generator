include("./generator_functions_DG.jl")

##### Set Spatial/Temporal Parameters Below

# Set time domain and discrete time spacing
Δt = 0.1
T_range = 0:Δt:1

# Create a grid and indexing for the spatial domain [xmin,xmax]x[ymin,ymax]
println("Setting up the grid...")
xmin, Δx, xmax = 0, 0.1, 3 
ymin, Δy, ymax = 0, 0.1, 2
d, grid = make_dict_grid(xmin, xmax, Δx, ymin, ymax, Δy)

##### Parameter Selection Complete

# Set the spatial diffusion parameter ϵ and an initial heuristic for the temporal diffusion strength a

ϵ, a_init = get_parameters(grid, T_range)

# Here is the value of a chosen to approximately match the leading spatial and temporal eigenvalues
a = 0.115

@time begin println("Making inflated generator...")
    # Create a vector of generators for each discrete time point
    Gvec = []
    for t ∈ T_range
        G = make_generator(d, grid, x -> dg_velocity(t, x), ϵ)
        push!(Gvec, G)
    end
    # Assemble individual generators into the inflated generator
    𝐆 = make_inflated_generator(Gvec, Δt, a)
end

println("Computing inflated generator eigenvalues...")
@time Λ, V = eigs(𝐆, which=:LR, nev=10, maxiter=100000)

# Plot of the spectrum
spectrumpicname = "./double gyre/Inflated Generator Eigenvalue Spectrum for the Double Gyre.png"
@time real_spat_inds = plot_spectrum_and_get_real_spatial_eigs(grid, Λ, V, spectrumpicname)

println("Plotting eigenvector time slices...")
# Plot slices of leading spatial eigenvector

index_to_plot = real_spat_inds[end]
time_slice_spacing = 1 # Plot every time_slice_spacing-th slice after start_date (time gap of time_slice_spacing*time_step)

picfilename = "./double gyre/Leading Spatial Eigenvector for the Double Gyre.png"
moviefilename = "./double gyre/Movie of leading spatial inflated generator eigenvector for the Double Gyre.gif"
@time plot_slices(real.(V), index_to_plot, time_slice_spacing, grid, T_range, :RdBu, picfilename, moviefilename)
# We have to use real() when producing plots for V or else an error will be thrown when attempting to plot the eigenvectors, even if imag(V[:,k]) = zeros(size(V,1)), as V is a matrix of complex type.

# Calculate SEBA Vectors from the leading two eigenvectors
println("Computing SEBA vectors...")
@time Σ, ℛ = SEBA(real.(V[:, real_spat_inds])) # Again, we must take the real part of V or an error will be thrown when running SEBA().
println("The respective SEBA vector minima are ", minimum(Σ, dims=1))

# Plot individual SEBA vector(s), followed by the maximum of the two
println("Plotting SEBA vector time slices...")

for index_to_plot = 1:length(real_spat_inds)
    picfilename = "./double gyre/SEBA vector $index_to_plot for the Double Gyre.png"
    moviefilename = "./double gyre/Movie of SEBA vector $index_to_plot for the Double Gyre.gif"
    @time plot_slices(Σ, index_to_plot, time_slice_spacing, grid, T_range, :Reds, picfilename, moviefilename)
end

println("Plotting time slices of SEBA vector maxima...")

Σ_max = maximum(Σ,dims=2)

index_to_plot = 1
picfilename = "./double gyre/Maxima of SEBA vectors for the Double Gyre.png"
moviefilename = "./double gyre/Movie of SEBA vector maxima for the Double Gyre.gif"
@time plot_slices(Σ_max, index_to_plot, time_slice_spacing, grid, T_range, :Reds, picfilename, moviefilename)

# Save the results to HDF5 and JLD2 files 
# Data to save: Vectors of grid ranges in x and y (or the entire grid struct in JLD2), time range vector, time slice spacing for plots, eigenvalues and eigenvectors of the inflated generator and SEBA vectors
println("Saving variables...")
filename = "./double gyre/InfGen_Results_SwitchingDoubleGyre"
@time save_results(grid, T_range, time_slice_spacing, Λ, V, Σ, filename)