using HDF5
include("generator_functions_DG.jl")

#WE DON'T NEED THIS BECAUSE IT IS IN GENERATOR_FUNCTIONS
include("plot_slices.jl")

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
# Value of F_median recorded: 0.7184901312589542

# Set the spatial diffusion parameter ϵ
ϵ = sqrt(0.1 * F_median * (grid.Δ_x))
println("The calculated ϵ value is... $ϵ")
# Value of ϵ recorded: 0.05360933244348242

# Set temporal diffusion parameter strength
a = sqrt(1.1 * F_median * (grid.Δ_x)) / 3
println("The initial a value is... $a")
# Value of a recorded: 0.05926734699215261
a = 0.1
# Value of a used: 0.1

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

#MISSING PLOT OF THE SPECTRUM

println("Plotting eigenvector time slices...")
#WHY ARE WE PLOTTING THE THIRD EIGENVECTOR AND NOT THE SECOND?
@time plot_slices(V, grid, 3)

# Calculate SEBA Vectors from the leading two eigenvectors
# UNCLEAR HOW VECTORS 1 AND 3 ARE THE LEADING TWO EIGENVECTORS
println("Computing SEBA vectors...")
seba_inds = [1, 3]
@time Σ, ℛ = SEBA(real.(V[:, seba_inds]))
println("The respective SEBA vector minima are ", minimum(Σ, dims=1))

# WE DON'T NEED SPECIAL SEBA PLOTTING CODE, JUST INPUT SEBA VECTORS INTO PLOT_SLICES
println("Plotting SEBA vector time slices...")
@time plot_SEBA(Σ, grid, 0) # For a max(SEBA) plot, insert 0 for vecnum

# Save the results to an HDF5 file 
println("Saving variables...")
#name_save_file = "InfGen_Results_DG_" * string(year(time_now)) * lpad(month(time_now), 2, "0") * lpad(day(time_now), 2, "0") * "_" * lpad(hour(time_now), 2, "0") * lpad(minute(time_now), 2, "0") * ".h5"
name_save_file = "InfGen_Results_SwitchingDoubleGyre.h5"
#A DESCRIPTION OF WHAT IS ACTUALLY BEING SAVED WOULD BE HELPFUL, SINCE IT SEEMS YOU DON'T SAVE THE INPUTS
@time save_results(grid, Λ, V, Σ, name_save_file)