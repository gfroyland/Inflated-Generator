using Statistics
using Dates
using HDF5
include("generator_functions_DG.jl")

# x and y arrays remain the same at each time step

Gvec = []

println("Setting up the grid...")
# 2/50
xmin, Δx, xmax = 0, 0.04, 3
ymin, Δy, ymax = 0, 0.04, 2

Δt = 0.05
T_range = 0:Δt:1

d, grid = make_dict_grid(xmin, xmax, Δx, ymin, ymax, Δy)

r(t) = (1 / 2) * (1 + tanh(10 * (t - (1 / 2))))
α(t) = (1 - 2 * r(t)) / (3 * (r(t) - 2) * (r(t) + 1))
β(t) = (2 - 9 * α(t)) / 3

F(t, x) = [(-π / 2) * sin(π * (α(t) * x[1]^2 + β(t) * x[1])) * cos(π * x[2] / 2), (2 * α(t) * x[1] + β(t)) * cos(π * (α(t) * x[1]^2 + β(t) * x[1])) * sin(π * x[2] / 2)]

F_median = median(norm(F(t, x)) for t ∈ T_range for x ∈ grid.centres)
println("The median of the speeds is... $F_median")
# 0.7184901312589542

ϵ = sqrt(0.1*F_median*(grid.Δ_x))
println("The calculated ϵ value is... $ϵ")
# 0.05360933244348242

@showprogress for t ∈ T_range

    #noise of 1 is reasonable since the integral is over side face of order 1
    #and the vector field norm ranges from 0 to 20.
    G = make_generator(d, grid, x -> F(t, x), ϵ)

    push!(Gvec, G)

end

Gᴰ = make_dynamic_generator(Gvec)

#compute 10 eigenvalues/vectors of largest real part
println("Solving eigenproblem...")
@time λ, v = eigs(Gᴰ, which=:LR, nev=10, maxiter=100000)
#for small problems, it may be better to convert G to a dense matrix and compute all eigenvalues/vectors
#@time λ, v = eigen(Matrix(G), sortby=t -> abs(t))

#calculate 6 SEBA vectors
numseba = 2
s, R = SEBA(real.(v[:, 1:numseba]))
println("Calculated $numseba SEBA vectors, respective minimum values are ", minimum(s, dims=1))

println("Plotting...")
@time plot_things(grid, λ, v, s)

a = sqrt(1.1*F_median*(grid.Δ_x))/3
println("The initial a value is... $a")
# 0.05926734699215261
# a = 0.125

println("Making inflated generator...")
@time 𝐆 = make_inflated_generator(Gvec, Δt, a)

println("Computing inflated eigenvalues...")
@time Λ, V = eigs(𝐆, which=:LR, nev=10, maxiter=100000)

println("Plotting slices...")
#@time plot_slices(V, grid, 3)
@time plot_spatemp_IDL(grid, Λ, V)
@time plot_9vecs_IDL(grid, Λ, V)

#=
for a ∈ 0.08:0.001:0.1

    println("Making inflated generator...")
    @time 𝐆 = make_inflated_generator(Gvec, Δt, a)

    println("Computing inflated eigenvalues...")
    @time Λ, V = eigs(𝐆, which=:LR, nev=10, maxiter=100000)

    println("Plotting slices...")
    @time plot_spatemp_IDL(grid, Λ, V)

end
=#

seba_inds = [1 ; 2]

Σ, ℛ = SEBA(real.(V[:, seba_inds]))
println("The respective SEBA vector minima are ", minimum(Σ, dims=1))

@time plot_SEBA_IDL(grid, Σ)

time_now = now()
name_save_file = "IDL_Results_DG_" * string(year(time_now)) * lpad(month(time_now), 2, "0") * lpad(day(time_now), 2, "0") * "_" * lpad(hour(time_now), 2, "0") * lpad(minute(time_now), 2, "0") * ".h5"
@time save_results(grid, λ, v, s, Λ, V, Σ, name_save_file)

#streamplot code...very slow with interpolated vector fields
#=
using CairoMakie
#V2(x::Point2f) = Point2f(F(x)[1], F(x)[2])
#V(x::Point2{T}) where T = Point2f(F(x)[1],F(x)[2])
VF(x::Point2{T}) where {T} = Point2f(F(x)[1], F(x)[2])
@time fig = streamplot(VF, lonmin .. lonmax, latmin .. latmax, colormap=:Blues)
#save("streamplot_1_august_2003.png", fig)

#arrows(grid.lonrange, grid.latrange, vec ∘ F)
=#


# 𝕄, 𝔸