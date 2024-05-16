using HDF5, JLD2
include("generator_functions_DG.jl")

# Set time domain and discrete time spacing
Δt = 0.05
T_range = 0:Δt:1

# Create a grid and indexing for the spatial domain [xmin,xmax]x[ymin,ymax]
println("Setting up the grid...")
xmin, Δx, xmax = 0, 0.1, 3 
ymin, Δy, ymax = 0, 0.1, 2
d, grid = make_dict_grid(xmin, xmax, Δx, ymin, ymax, Δy)

# Define the time-dependent switching double gyre vector field F(t,x) see [Atnip/Froyland/Koltai, 2024]
r(t) = (1 / 2) * (1 + tanh(10 * (t - (1 / 2))))
α(t) = (1 - 2 * r(t)) / (3 * (r(t) - 2) * (r(t) + 1))
β(t) = (2 - 9 * α(t)) / 3
F(t, x) = [(-π / 2) * sin(π * (α(t) * x[1]^2 + β(t) * x[1])) * cos(π * x[2] / 2), (2 * α(t) * x[1] + β(t)) * cos(π * (α(t) * x[1]^2 + β(t) * x[1])) * sin(π * x[2] / 2)]

F_median = median(norm(F(t, x)) for t ∈ T_range for x ∈ grid.centres)
println("The median of the speeds is... $F_median")

# Set the spatial diffusion parameter ϵ
ϵ = sqrt(0.1 * F_median * (grid.Δ_x))
println("The calculated ϵ value is... $ϵ")

# Set temporal diffusion parameter strength
L_max_x = (grid.x_max - grid.x_min)
L_max_y = (grid.y_max - grid.y_min)
a = sqrt(1.1 * F_median * (grid.Δ_x)) / (max(L_max_x,L_max_y)) # Initial heuristic for a
println("The initial a value is... $a")
a = 0.115 # a chosen to approximately match the leading eigenvalues

@time begin println("Making inflated generator...")
    # Create a vector of generators for each discrete time point
    Gvec = []
    for t ∈ T_range
        G = make_generator(d, grid, x -> F(t, x), ϵ)
        push!(Gvec, G)
    end
    # Assemble individual generators into the inflated generator
    𝐆 = make_inflated_generator(Gvec, Δt, a)
end

println("Computing inflated generator eigenvalues...")
@time Λ, V = eigs(𝐆, which=:LR, nev=10, maxiter=100000)

# PLOT OF THE SPECTRUM
@time plot_spectrum(grid, Λ, V)

println("Plotting eigenvector time slices...")
# Plot slices of leading spatial eigenvector (V_2)
# Make sure to ℓ^2-normalise V before sending it through

V_norm = stack(normalize.(eachcol(V))) * sqrt(size(V, 1))
vecnum = 2
@time plot_slices(real.(V_norm), vecnum, grid, T_range, :RdBu) 
# We have to use real() when producing plots for V, or else an error will be thrown when attempting to plot the eigenvectors, even if imag(V[:,vecnum]) = zeros(size(V,1)), as V is a matrix of complex type.

# Calculate SEBA Vectors from the leading two eigenvectors
println("Computing SEBA vectors...")
seba_inds = [1, 2]
@time Σ, ℛ = SEBA(real.(V[:, seba_inds]))
println("The respective SEBA vector minima are ", minimum(Σ, dims=1))

# Plot individual SEBA vectors, followed by the maximum of the two
println("Plotting SEBA vector time slices...")

sebanum = 1
@time plot_slices(Σ, sebanum, grid, T_range, :Reds)

Σ_max = maximum(Σ,dims=2)
@time plot_slices(Σ_max, 1, grid, T_range, :Reds)

# Save the results to HDF5 and JLD2 files 
# Data to save: Vectors of grid ranges in x and y, time range vector, eigenvalues and eigenvectors of the inflated generator and SEBA vectors
println("Saving variables...")
name_save_file = "InfGen_Results_SwitchingDoubleGyre"
@time save_results(grid, T_range, Λ, V, Σ, name_save_file)