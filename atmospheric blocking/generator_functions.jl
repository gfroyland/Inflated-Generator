using Dates, Plots, LinearAlgebra, Interpolations, QuadGK, SparseArrays, Arpack, Statistics, ProgressMeter, HDF5, JLD2

#create a data structure for the grid;  it might be expanded later to include e.g. G, spectrum, evecs, etc... or other items
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

include("../SEBA.jl")

#create a dictionary to do the indexing we want and a grid struct
function make_dict_grid(lonmin, lonmax, lonspacing, latmin, latmax, latspacing)
    lonrange = lonmin+lonspacing/2:lonspacing:lonmax-lonspacing/2
    latrange = latmin+latspacing/2:latspacing:latmax-lonspacing/2
    lonrange = round.(lonrange, digits=6)
    latrange = round.(latrange, digits=6)
    i = 0
    temparray = []
    for lon ∈ lonrange
        for lat ∈ latrange
            i += 1
            push!(temparray, ([lon, lat], i))
        end
    end
    d = Dict(temparray)
    centres = Tuple.(collect(keys(d)))
    grid = Grid(centres, lonrange, latrange, lonmin, lonmax, latmin, latmax, lonspacing, latspacing)
    return d, grid
end

function get_linear_interpolant(lons_data, lats_data, u_data, v_data)

    # Here, lons_data represents longitude and lats_data represents latitude
    # There should be no need to thin the data from here on in

    # The reason that reverse() is used below is because the latitudinal coordinate vector
    # supplied by ERA5 provides the latitudinal range in descending rather than ascending order.
    # This vector (and the u/v component matrices in the second dimension) therefore need to be reversed, or else the convenience
    # constructors used below will (in the same fashion as MATLAB's griddedInterpolant) throw an error 
    # because the grid coordinates have not been supplied in ascending order.

    Iplt_z = linear_interpolation((lons_data, reverse(lats_data)), reverse(u_data, dims=2))
    Iplt_m = linear_interpolation((lons_data, reverse(lats_data)), reverse(v_data, dims=2))

    return Iplt_z, Iplt_m
end

function make_generator(d, grid, F, ϵ)

    # The strategy is to follow Froyland/Junge/Koltai (2013).  
    # Spatial coordinates are in degrees lon/lat (info contained in d and grid)
    # Velocities (the output of F) are in m.s⁻¹
    # ϵ is a diffusion parameter, making a crude finite-difference diffusion on a 5-point stencil 
    
    #conversion factor from degrees to metres:  mean circumference of Earth at the equator is 40075000m and there are 2π radians (360 degrees) in a circle, so to obtain metres from degrees at the equator we multiply by 40075000/360
    deg2metr = 40075000/360

    #create list of box centres
    centres = collect(keys(d))
    
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
    tol = 1e-2  #hard coded for now, but ultimately should be adapted to entry sizes of 
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

function make_dynamic_generator(Gvec)

    Gᴰ = mean(Gvec)
    return Gᴰ

end

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

function plot_spectrum(grid, Λ, V)

    # Calculate Means of Eigenvector Variances 

    N = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / N)

    K = size(V, 2)

    meanvariance = [mean([var(V[(t-1)*N+1:t*N, k]) for t = 1:T]) for k = 1:K]

    # Plot the Spectrum, distinguishing spatial eigenvalues from temporal ones
    
    spat_inds = findall(x->x>1e-10,meanvariance)
    temp_inds = findall(x->x<1e-10,meanvariance)

    # Trivial Λ_1 should be plotted as a spatial eigenvalue, but meanvariance[1] ≈ 0, alleviorate this before plotting

    popfirst!(temp_inds) 
    append!(spat_inds,1)

    scatter(Λ[spat_inds], label="Spatial Λ_k", shape=:circle, mc=:blue, title="$(length(Λ)) eigenvalues with largest real part, a = $a", xlabel="Re(Λ_k)", ylabel="Im(Λ_k)")
    scatter!(Λ[temp_inds], label="Temporal Λ_k", shape=:xcross, mc=:red, msw=4)
    xlabel!("Re(Λ_k)")
    display(ylabel!("Im(Λ_k)"))

end

function plot_slices(V, vecnum, time_step, grid, date_range, col_scheme, figlayout)

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = length(date_range)

    #create a T-vector of time-slices (copies of space)
    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    # find a common colour range
    col_lims = (minimum((V[:, vecnum])),maximum((V[:, vecnum])))

    # create an animation of frames of the eigenvector
    anim = @animate for t = 1:time_step:T
        title_now = Dates.format(date_range[t], "dd/mm HH:MM")
        contourf(grid.lonrange, grid.latrange, reshape(sliceV[t][:, vecnum], length(grid.latrange), length(grid.lonrange)), clims=col_lims, c=col_scheme, xlabel="̊ E", ylabel="̊ N", title=title_now, linewidth=0, levels=100)
    end
    display(gif(anim, fps=10))

    # plot individual time frames
    fig = []
    for t = 1:time_step:T
        title_now = Dates.format(date_range[t], "dd/mm HH:MM")
        push!(fig, contourf(grid.lonrange, grid.latrange, reshape(sliceV[t][:, vecnum], length(grid.latrange), length(grid.lonrange)), clims=col_lims, c=col_scheme, title=title_now, linewidth=0, levels=100, aspectratio=1, legend=:none))
    end
    display(plot(fig..., layout=figlayout))

end

function save_results(grid, date_range, Λ, V, Σ, filename)
    
    # Save data to JLD2 file
    filename_JLD2 = filename * ".jld2"
    jldsave(filename_JLD2; grid, date_range, Λ, V, Σ)

    # Save data to HDF5 file
    filename_HDF5 = filename * ".h5"
    file_ID = h5open(filename_HDF5, "w")

    file_ID["lonrange"] = grid.lonrange
    file_ID["latrange"] = grid.latrange

    date_range_str = String[]
    for i = 1:length(date_range)
        push!(date_range_str,Dates.format(date_range[i],"dd/mm/yy HH:MM"))
    end

    file_ID["date_range"] = collect(date_range_str) # The collect() function must be used or an error will be thrown

    file_ID["Eigvals_Real"] = real.(Λ)
    file_ID["Eigvals_Imag"] = imag.(Λ)

    file_ID["Eigvecs_Real"] = real.(V)
    file_ID["Eigvecs_Imag"] = imag.(V)

    file_ID["SEBA"] = Σ

    close(file_ID)

end
