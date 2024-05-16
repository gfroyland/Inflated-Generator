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

include("SEBA.jl")

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
    jldsave(filename_JLD2; grid.lonrange, grid.latrange, date_range, Λ, V, Σ)

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
#=
function quiverplus(xrange, yrange, F, sparsity)
    # e.g. F(x,y)=[-y, x]
    # quiverplus(-1:0.1:1, -1:0.1:1, F, 6)
    # sparsity of 6 means keeping each 6th point in the 2D grid (to make the plot less busy)

    vF(x) = F(x...)
    centres = [[x, y] for x ∈ xrange for y ∈ yrange]
    sparsecentres = sort(centres)[1:sparsity:end]
    sparsecentremat = hcat(sparsecentres...)'
    meannorm = mean(norm.(vF.(centres)))  #mean vector norm
    scale = norm([sort(xrange)[2] - sort(xrange)[1], sort(yrange)[2] - sort(yrange)[1]]) * (√sparsity) / meannorm #scale quiver vectors by multiplicative factor "scale", computed so that the average-length vector doesn't overlap
    Fs(x, y) = scale * F(x, y)
    display(quiver(sparsecentremat[:, 1], sparsecentremat[:, 2], quiver=Fs, arrow=:filled, xlabel="x", ylabel="y"))
end

function plot_vector_field(d, grid, F)

    #plots windspeed field
    #FF(x, y) = F([x, y])
    #display(surface(grid.lonrange, grid.latrange, norm ∘ FF, xlabel="lon", ylabel="lat", zlabel="wind speed (units?)", c=:jet))
    #display(Plots.contourf(grid.lonrange, grid.latrange, norm ∘ FF, xlabel="lon", ylabel="lat", title="wind speed (units)?", aspectratio=1, c=:jet))

    #plot quiver of windspeed field, with scaling of arrows inspired by the matlab autoscaling
    centres = collect(keys(d))
    sparsity = 7  #include only every sparsity^{th} grid point, ideally sparsity is coprime with the number of grid points in either lon or lat direction
    sparsecentres = sort(centres)[1:sparsity:end]
    sparsecentremat = hcat(sparsecentres...)'
    #vf=F.(sparsecentres)
    #vfmat=vcat(vf...)
    meannorm = mean(norm.(F.(sparsecentres)))  #mean vector norm
    scale = norm([grid.lonspacing, grid.latspacing]) * (√sparsity) / meannorm #scale quiver vectors by multiplicative factor "scale", computed so that the average-length vector doesn't overlap
    vF(x, y) = scale * FF(x, y)
    display(quiver(sparsecentremat[:, 1], sparsecentremat[:, 2], quiver=vF, title="wind field", arrow=:filled, aspectratio=1, xlabel="lat", ylabel="lon"))

    #plot a coloured spy plot of G
    #plotlyjs()
    #display(heatmap(G, c=:RdBu, yflip=true, title="heatmap of matrix G"))
    #gr()

end

function plot_things(grid, λ, v, s)

    #plot spectrum
    display(Plots.scatter(λ, title="$(length(λ)) eigenvalues with largest real part"))

    #plots image of an eigenvector in lon,lat coordinates
    #display(heatmap(grid.lonrange, grid.latrange, real.(reshape(v[:, evecnum], length(grid.latrange), length(grid.lonrange))), xlabel="longitude", ylabel="latitude", title="2D image of eigenvector number $evecnum", c=:RdBu, aspectratio=1))

    #plots leading 9 eigenvectors
    # aspectratio = 1
    p = []
    for k = 1:9
        push!(p, heatmap(grid.lonrange, grid.latrange, real.(reshape(v[:, k], length(grid.latrange), length(grid.lonrange))), c=:RdBu, title="eigenvector $k"))
    end
    display(Plots.plot(p..., layout=(3, 3)))

    #plots leading seba vectors
    pseba = []
    for k = 1:size(s)[2]
        push!(pseba, heatmap(grid.lonrange, grid.latrange, real.(reshape(s[:, k], length(grid.latrange), length(grid.lonrange))), c=:Reds, title="SEBA vector $k"))
    end
    display(Plots.plot(pseba..., layout=(3, 3)))

    #plot max of seba vectors
    smax = maximum(s, dims=2)
    display(heatmap(grid.lonrange, grid.latrange, real.(reshape(smax, length(grid.latrange), length(grid.lonrange))), title="superposition of $(size(s)[2]) SEBA vector(s)", c=:Reds, xlabel="lon", ylabel="lat"))

    #imshow(real.(reshape(smax, length(grid.latrange), length(grid.lonrange))), grid.lonmin:grid.lonspacing:grid.lonmax,  grid.latmin:grid.latspacing:grid.latmax, projection=:Robinson, coast=true, colorbar=true, cmap=:red2green, fmt=:png, title="superposition of $(size(s)[2]) SEBA vectors")
end

function plot_things_IDL(grid, Λ, V, numseba)

    #plot inflated spectrum
    display(Plots.scatter(Λ, title="$(length(Λ)) eigenvalues with largest real part, a = $a"))

    # Prepare the slices of nine eigenvectors at each time step 

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)
    xy = [[x, y] for t = 1:T for x ∈ grid.lonrange for y ∈ grid.latrange]
    xy = hcat(xy...)'

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    𝒫_Σ = [] # For the superposition of SEBA Vectors

    for τ ∈ 1:T

        #plots leading 9 eigenvectors
        P = []
        for k = 1:9
            push!(P, heatmap(grid.lonrange, grid.latrange, real.(reshape(sliceV[τ][:, k], length(grid.latrange), length(grid.lonrange))), c=:RdBu, title="v_$k, D$(div(τ-1,4)), $(((τ-1)%4)*6)h", aspectratio=1))
        end
        display(Plots.plot(P..., layout=(3, 3)))

        #calculate and plot leading seba vectors
        if numseba == 1
            Σ, R = SEBA(real.(sliceV[τ][:, 9]))
        else
            Σ, R = SEBA(real.(sliceV[τ][:, 1:numseba]))
        end
        P_Σ = []
        for k = 1:size(Σ)[2]
            push!(P_Σ, heatmap(grid.lonrange, grid.latrange, real.(reshape(Σ[:, k], length(grid.latrange), length(grid.lonrange))), c=:Reds, title="SEBA($k), D$(div(τ-1,4)), $(((τ-1)%4)*6)h", aspectratio=1))
        end
        display(Plots.plot(P_Σ..., layout=(3, 3)))

        #plot max of seba vectors
        Σ_max = maximum(Σ, dims=2)
        push!(𝒫_Σ, heatmap(grid.lonrange, grid.latrange, real.(reshape(Σ_max, length(grid.latrange), length(grid.lonrange))), c=:Reds, title="$numseba xSEBA, D$(div(τ-1,4)), $(((τ-1)%4)*6)h", aspectratio=1))

        #imshow(real.(reshape(smax, length(grid.latrange), length(grid.lonrange))), grid.lonmin:grid.lonspacing:grid.lonmax,  grid.latmin:grid.latspacing:grid.latmax, projection=:Robinson, coast=true, colorbar=true, cmap=:red2green, fmt=:png, title="superposition of $(size(s)[2]) SEBA vectors")

    end

    # Be advised, the layout dimensions below will have to change with T 
    # This will have to be done manually for now...
    #display(Plots.plot(𝒫_Σ..., layout=(7, 3)))

end

function plot_eigvecs_IDL(grid, Λ, V, vecnum)

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    cmax = maximum(abs.(real(V[:, vecnum])))

    anim = @animate for τ ∈ 1:T

        display(Plots.contourf(grid.lonrange, grid.latrange, reshape(real.(sliceV[τ][:, vecnum]), length(grid.latrange), length(grid.lonrange)), clims=(-cmax, cmax), c=:RdBu, xlabel="x", ylabel="y", title="v_$vecnum, Day $(div(τ-1,4)), $(((τ-1)%4)*6):00", linewidth=0, levels=100))

    end

    gif(anim, "EuroBlock_3DayExt_Vec8_a3p5_24Jan24.gif", fps=2)
    gif(anim, "EuroBlock_3DayExt_Vec8_a3p5_24Jan24.mp4", fps=2)

end

function plot_9vecs_IDL(grid, V, date_range)

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    cmax = maximum(real(V), dims=1)
    cmin = minimum(real(V), dims=1)
    #cmax[6] = 0.018 # Use this for East Block, No Extension
    #cmax[7] = 0.0001 # Use this for East Block, No Extension
    
    for τ ∈ 1:T
    #anim = @animate for τ ∈ 1:T

        P = []
        for k = 1:9
            title_now = "v_$k, " * Dates.format(date_range[τ], "dd/mm HH:MM")
            push!(P, Plots.contourf(grid.lonrange, grid.latrange, reshape(real.(sliceV[τ][:, k]), length(grid.latrange), length(grid.lonrange)), clims=(cmin[k], cmax[k]), c=:RdBu, xlabel="x", ylabel="y", title=title_now, linewidth=0, levels=100))
        end
        display(Plots.plot(P..., layout=(3, 3)))
    end

    #gif(anim, "EuroBlock_3DayExt_9EigVecs_a3p3_05Feb24.gif", fps=2)
    #gif(anim, "EuroBlock_East_3DayExtension_9EigVecs_a0p7_16Feb24.mp4", fps=2)

end

function plot_Nvecs_IDL(grid, V, date_range)

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    cmax = maximum(real(V), dims=1)
    cmin = minimum(real(V), dims=1)
    #cmax[6] = 0.018 # Use this for East Block, No Extension
    #cmax[7] = 0.0001 # Use this for East Block, No Extension

    N = ceil(Int64,size(V,2)/9)

    for n ∈ 1:N

        ind_start = (9*(n-1))+1
        ind_end = 9*n

        if (ind_end > size(V,2))
            ind_end = size(V,2)
        end
    
    for τ ∈ 1:T
    #anim = @animate for τ ∈ 1:T

        P = []
        for k = ind_start:ind_end
            title_now = "v_$k, " * Dates.format(date_range[τ], "dd/mm HH:MM")
            push!(P, Plots.contourf(grid.lonrange, grid.latrange, reshape(real.(sliceV[τ][:, k]), length(grid.latrange), length(grid.lonrange)), clims=(cmin[k], cmax[k]), c=:RdBu, xlabel="x", ylabel="y", title=title_now, linewidth=0, levels=100))
        end
        display(Plots.plot(P..., layout=(3, 3)))
    end

    #vidname = "EuroBlock_East_3DayExtension_EigVecs_" * lpad(ind_start, 1, "0") * "to" * lpad(ind_end, 1, "0") * "_a0p7_12Mar24.mp4"
    
    #println(vidname)
    #gif(anim, vidname, fps=2)

    end

end

function plot_specvecs_IDL(grid, V, date_range, eigs_to_plot)

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    cmax = maximum(real(V), dims=1)
    cmin = minimum(real(V), dims=1)
    #cmax[6] = 0.018 # Use this for East Block, No Extension
    #cmax[7] = 0.0001 # Use this for East Block, No Extension

    #for τ ∈ 1:T
    anim = @animate for τ ∈ 1:T

        P = []
        for k = 1:length(eigs_to_plot) # eigs_to_plot must be of length 9 or less, and must not contain indices larger than size(V,2)
            vecnow = eigs_to_plot[k]
            title_now = "v_$vecnow, " * Dates.format(date_range[τ], "dd/mm HH:MM")
            push!(P, Plots.contourf(grid.lonrange, grid.latrange, reshape(real.(sliceV[τ][:, vecnow]), length(grid.latrange), length(grid.lonrange)), clims=(cmin[vecnow], cmax[vecnow]), c=:RdBu, xlabel="x", ylabel="y", title=title_now, linewidth=0, levels=100))
        end
        display(Plots.plot(P..., layout=(3, 3)))
    end

    vidname = "EuroBlock_West_5DayExtension_18Eigvecs_a1p2_11Mar24.mp4"
    
    println(vidname)
    gif(anim, vidname, fps=2)

end

function plot_9leftvecs_IDL(grid, V, date_range)

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    cmax = maximum(real(V), dims=1)
    cmin = minimum(real(V), dims=1)
    #cmax[6] = 0.018 # Use this for East Block, No Extension
    #cmax[7] = 0.0001 # Use this for East Block, No Extension
    # 0.008 (emergency scaling)
    
    for τ ∈ 1:T
    #anim = @animate for τ ∈ 1:T

        P = []
        for k = 1:9
            title_now = "v_$k, " * Dates.format(date_range[τ], "dd/mm HH:MM")
            push!(P, Plots.contourf(grid.lonrange, grid.latrange, reshape(real.(sliceV[τ][:, k]), length(grid.latrange), length(grid.lonrange)), clims=(cmin[k], cmax[k]), c=:RdBu, xlabel="x", ylabel="y", title=title_now, linewidth=0, levels=100))
        end
        display(Plots.plot(P..., layout=(3, 3)))
    end

    #gif(anim, "EuroBlock_3DayExt_9EigVecs_a3p3_05Feb24.gif", fps=2)
    #gif(anim, "EuroBlock_East_3DayExtension_9EigVecs_Transpose_a0p7_New_16Feb24.mp4", fps=2)

end

function plot_9leftvecs_log_IDL(grid, V, date_range)

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)

    logV = real.(V)

    Vpos = findall(x->x>0,logV)
    Vneg = findall(x->x<0,logV)

    logV[Vpos] = -log.(logV[Vpos])
    logV[Vneg] = log.(abs.(logV[Vneg]))

    sliceV = [logV[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    cmax = maximum(logV, dims=1)
    cmin = minimum(logV, dims=1)
    #cmax[6] = 0.018 # Use this for East Block, No Extension
    #cmax[7] = 0.0001 # Use this for East Block, No Extension
    
    #for τ ∈ 1:T
    anim = @animate for τ ∈ 1:T

        P = []
        for k = 1:9
            title_now = "u_$k, " * Dates.format(date_range[τ], "dd/mm HH:MM")
            push!(P, Plots.contourf(grid.lonrange, grid.latrange, reshape(sliceV[τ][:, k], length(grid.latrange), length(grid.lonrange)), clims=(cmin[k], cmax[k]), c=:RdBu, xlabel="x", ylabel="y", title=title_now, linewidth=0, levels=100))
        end
        display(Plots.plot(P..., layout=(3, 3)))
    end

    #gif(anim, "EuroBlock_3DayExt_9EigVecs_a3p3_05Feb24.gif", fps=2)
    gif(anim, "EuroBlock_East_3DayExtension_9EigVecs_Transpose_a0p7_Log_16Feb24.mp4", fps=2)

end

function plot_9morevecs_IDL(grid, V, date_range)

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    cmax = maximum(real(V), dims=1)
    cmin = minimum(real(V), dims=1)
    #cmax[6] = 0.018 # Use this for East Block, No Extension
    #cmax[7] = 0.0001 # Use this for East Block, No Extension
    
    for τ ∈ 1:T
    #anim = @animate for τ ∈ 1:T

        P = []
        for k = 10:18
            title_now = "v_$k, " * Dates.format(date_range[τ], "dd/mm HH:MM")
            push!(P, Plots.contourf(grid.lonrange, grid.latrange, reshape(real.(sliceV[τ][:, k]), length(grid.latrange), length(grid.lonrange)), clims=(cmin[k], cmax[k]), c=:RdBu, xlabel="x", ylabel="y", title=title_now, linewidth=0, levels=100))
        end
        display(Plots.plot(P..., layout=(3, 3)))
    end

    #gif(anim, "EuroBlock_3DayExt_9EigVecs_a3p3_05Feb24.gif", fps=2)
    #gif(anim, "EuroBlock_East_3DayExtension_9MoreEigVecs_a0p7_16Feb24.mp4", fps=2)

end

function make_SEBA_IDL(grid, V, seba_inds, date_range)

    # Prepare the slices of nine eigenvectors at each time step 
    numseba = length(seba_inds)
    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    Vecs_for_SEBA = zeros(size(V)[1], (2 * numseba))

    for η ∈ 1:numseba

        Vecs_for_SEBA[:, (2*η-1)] = real.(V[:, seba_inds[η]])

    end

    for τ ∈ 1:T

        for κ ∈ 1:numseba
            Vecs_for_SEBA[(τ-1)*spacelength.+(1:spacelength), (2*κ)] = norm(sliceV[τ][:, seba_inds[κ]]) .* ones(spacelength)
        end

    end

    Σ, ℛ = SEBA(Vecs_for_SEBA)
    println("The respective SEBA vector minima are ", minimum(Σ, dims=1))

    # Plot All SEBA Evolutions Once Finished

    
    sliceΣ = [Σ[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    for σ ∈ 1:size(Σ)[2]
        vidname = "EuroBlock_East_3DayExtension_Aug_NoCplx_SEBA" * lpad(σ, 1, "0") * "_a0p7_16Feb24.mp4"
        for t ∈ 1:T
        #anim = @animate for t ∈ 1:T

            title_now = "SEBA $σ, " * Dates.format(date_range[t], "dd/mm/yy HH:MM")
            display(Plots.contourf(grid.lonrange, grid.latrange, reshape(sliceΣ[t][:, σ], length(grid.latrange), length(grid.lonrange)), clims=(0, 1), c=:Reds, xlabel="x", ylabel="y", title=title_now, linewidth=0, levels=100))

        end
        #gif(anim, vidname, fps=2)
        println("SEBA Video $σ Complete.")
    end

    # Plot the Maxima of each SEBA Vector

    vidname = "EuroBlock_East_3DayExtension_Aug_NoCplx_SEBAMax_a0p7_16Feb24.mp4"
    for 𝒯 ∈ 1:T
    #anim = @animate for 𝒯 ∈ 1:T

        title_now = "Max SEBA, " * Dates.format(date_range[𝒯], "dd/mm/yy HH:MM")
        display(Plots.contourf(grid.lonrange, grid.latrange, reshape(maximum(sliceΣ[𝒯][:, :], dims=2), length(grid.latrange), length(grid.lonrange)), clims=(0, 1), c=:Reds, xlabel="x", ylabel="y", title=title_now, linewidth=0, levels=100))

    end
    #gif(anim, vidname, fps=2)
    println("Max SEBA Video Complete.")
    

    #return Σ
end

function plot_SEBA_IDL(grid, Σ, date_range)

    # Prepare the slices of nine eigenvectors at each time step 
    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(Σ)[1] / spacelength)

    sliceΣ = [Σ[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    # Plot All SEBA Evolutions

    for σ ∈ 1:size(Σ)[2]
        #vidname = "EuroBlock_East_3DayExtension_Upto13_SEBA" * lpad(σ, 1, "0") * "_a0p7_13Mar24.mp4"
        for t ∈ 1:T
        #anim = @animate for t ∈ 1:T

            title_now = "SEBA $σ, " * Dates.format(date_range[t], "dd/mm/yy HH:MM")
            display(Plots.contourf(grid.lonrange, grid.latrange, reshape(sliceΣ[t][:, σ], length(grid.latrange), length(grid.lonrange)), clims=(0, 1), c=:Reds, xlabel="x", ylabel="y", title=title_now, linewidth=0, levels=100))

        end
        #gif(anim, vidname, fps=2)
        println("SEBA Video $σ Complete.")
    end

    # Plot the Maxima of each SEBA Vector

    #vidname = "EuroBlock_East_3DayExtension_Upto13_SEBAMax_a0p7_13Mar24.mp4"
    for 𝒯 ∈ 1:T
    #anim = @animate for 𝒯 ∈ 1:T

        title_now = "Max SEBA, " * Dates.format(date_range[𝒯], "dd/mm/yy HH:MM")
        display(Plots.contourf(grid.lonrange, grid.latrange, reshape(maximum(sliceΣ[𝒯][:, :], dims=2), length(grid.latrange), length(grid.lonrange)), clims=(0, 1), c=:Reds, xlabel="x", ylabel="y", title=title_now, linewidth=0, levels=100))

    end
    #gif(anim, vidname, fps=2)
    println("Max SEBA Video Complete.")

end

function save_seba(grid, V, numseba, filename)

    file_ID = h5open(filename, "w")

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)
    xy = [[x, y] for t = 1:T for x ∈ grid.lonrange for y ∈ grid.latrange]
    xy = hcat(xy...)'

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    file_ID["longitude"] = grid.lonrange
    file_ID["latitude"] = grid.latrange

    for τ ∈ 1:T

        #calculate and plot leading seba vectors
        if numseba == 1
            Σ, R = SEBA(real.(sliceV[τ][:, 9]))
        else
            Σ, R = SEBA(real.(sliceV[τ][:, 1:numseba]))
        end

        Σ_max = maximum(Σ, dims=2)

        varname = "SEBA_" * lpad(τ, 2, "0")
        file_ID[varname] = real.(reshape(Σ_max, length(grid.latrange), length(grid.lonrange)))
    end

    close(file_ID)

end

function save_results(grid, λ, v, s, Λ, V, Σ, filename)

    file_ID = h5open(filename, "w")

    file_ID["longitude"] = grid.lonrange
    file_ID["latitude"] = grid.latrange

    s_max = maximum(s, dims=2)
    Σ_max = maximum(Σ, dims=2)

    file_ID["DL_Lambda"] = real.(λ)
    for β ∈ 1:size(v, 2)
        varname = "DL_v_" * lpad(β, 2, "0")
        file_ID[varname] = real.(reshape(v[:, β], length(grid.latrange), length(grid.lonrange)))
    end
    for γ ∈ 1:size(s,2)
        varname = "DL_SEBA_" * lpad(γ, 2, "0")
        file_ID[varname] = real.(reshape(s[:, γ], length(grid.latrange), length(grid.lonrange)))
    end
    file_ID["DL_SEBA_Max"] = real.(reshape(s_max, length(grid.latrange), length(grid.lonrange)))

    file_ID["IDL_Lambda"] = real.(Λ)

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    sliceΣ = [Σ[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    for τ ∈ 1:T

        for φ ∈ 1:size(V, 2)
            varname = "IDL_v_" * lpad(φ, 2, "0") * "_t" * lpad(τ, 2, "0")
            file_ID[varname] = real.(reshape(sliceV[τ][:, φ], length(grid.latrange), length(grid.lonrange)))
        end
        for ψ ∈ 1:size(Σ, 2)
            varname = "IDL_SEBA_" * lpad(ψ, 2, "0") * "_t" * lpad(τ, 2, "0")
            file_ID[varname] = real.(reshape(sliceΣ[τ][:, ψ], length(grid.latrange), length(grid.lonrange)))
        end

        varname = "IDL_SEBA_Max_t" * lpad(τ, 2, "0")
        file_ID[varname] = real.(reshape(Σ_max[(τ-1)*spacelength.+(1:spacelength)], length(grid.latrange), length(grid.lonrange)))
    end

    close(file_ID)

end

function save_results_DL(grid, λ, v, s, filename)

    file_ID = h5open(filename, "w")

    file_ID["longitude"] = grid.lonrange
    file_ID["latitude"] = grid.latrange

    s_max = maximum(s, dims=2)

    file_ID["DL_Lambda_Real"] = real.(λ)
    file_ID["DL_Lambda_Imag"] = imag.(λ)
    for β ∈ 1:size(v, 2)
        varname = "DL_v_" * lpad(β, 2, "0")
        file_ID[varname] = real.(reshape(v[:, β], length(grid.latrange), length(grid.lonrange)))
    end
    for γ ∈ 1:size(s,2)
        varname = "DL_SEBA_" * lpad(γ, 2, "0")
        file_ID[varname] = real.(reshape(s[:, γ], length(grid.latrange), length(grid.lonrange)))
    end
    file_ID["DL_SEBA_Max"] = real.(reshape(s_max, length(grid.latrange), length(grid.lonrange)))

    close(file_ID)

end

function save_results_IDL(grid, Λ, V, Σ, filename)

    file_ID = h5open(filename, "w")

    file_ID["longitude"] = grid.lonrange
    file_ID["latitude"] = grid.latrange

    Σ_max = maximum(Σ, dims=2)

    file_ID["IDL_Lambda_Real"] = real.(Λ)
    file_ID["IDL_Lambda_Imag"] = imag.(Λ)

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    sliceΣ = [Σ[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    for τ ∈ 1:T

        for φ ∈ 1:size(V, 2)
            varname = "IDL_v_" * lpad(φ, 2, "0") * "_t" * lpad(τ, 2, "0")
            file_ID[varname] = real.(reshape(sliceV[τ][:, φ], length(grid.latrange), length(grid.lonrange)))
        end
        for ψ ∈ 1:size(Σ, 2)
            varname = "IDL_SEBA_" * lpad(ψ, 2, "0") * "_t" * lpad(τ, 2, "0")
            file_ID[varname] = real.(reshape(sliceΣ[τ][:, ψ], length(grid.latrange), length(grid.lonrange)))
        end

        varname = "IDL_SEBA_Max_t" * lpad(τ, 2, "0")
        file_ID[varname] = real.(reshape(Σ_max[(τ-1)*spacelength.+(1:spacelength)], length(grid.latrange), length(grid.lonrange)))
    end

    close(file_ID)

end

function plot_slices(Λ, V, grid, vecnum, a)

    #X=[(t,x,y) for t=1:T  for x∈grid.lonrange for y∈grid.latrange ]
    #scatter(X,zcolor=real(V[:,2]))

    #plot inflated spectrum
    display(Plots.scatter(Λ, title="$(length(Λ)) eigenvalues with largest real part, a = $a"))

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)
    xy = [[x, y] for t = 1:T for x ∈ grid.lonrange for y ∈ grid.latrange]
    xy = hcat(xy...)'

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]
    slicenorm = [norm(sliceV[t][:, k]) for t = 1:T, k = 1:min(5, size(V)[2])]
    slicemean = [mean(sliceV[t][:, k]) for t = 1:T, k = 1:min(5, size(V)[2])]
    display(Plots.plot(Plots.plot(real(slicenorm), title="a = $a"), Plots.plot(real(slicemean)), layout=(2, 1)))

    slicenorm_2 = [norm(sliceV[t][:, k]) for t = 1:T, k = 6:10]
    slicemean_2 = [mean(sliceV[t][:, k]) for t = 1:T, k = 6:10]
    display(Plots.plot(Plots.plot(real(slicenorm_2), title="a = $a", labels=["y6" "y7" "y8" "y9" "y10"]), Plots.plot(real(slicemean_2), labels=["y6" "y7" "y8" "y9" "y10"]), layout=(2, 1)))

    #anim = @animate for t = 1:T
    #    display(heatmap(grid.lonrange, grid.latrange, reshape(real.(sliceV[t][:, vecnum]), length(grid.latrange), length(grid.lonrange)), c=:RdBu, xlabel="lon", ylabel="lat", title="Day $(div(t-1,2)), $(((t-1)%2)*12):00"))
    #display(imshow(reshape(real.(sliceV[t][:, vecnum]), length(grid.latrange), length(grid.lonrange)), grid.lonmin:grid.lonspacing:grid.lonmax,  grid.latmin:grid.latspacing:grid.latmax, projection=:Robinson, coast=true, colorbar=true, cmap=:red2green, fmt=:png, title="Day $t"))

    #imshow(-G, -5:0.5:42, 45:0.5:75, projection=:Robinson, coast=true, colorbar=true, fmt=:png, cmap=:red2green)
    #also cmap: panoply, hot, bam, vik, balance, buda, haline, matter, amp, oslo, gebco, curl
    #end

    #   gif(anim, "ecmwf-24jul-4aug-daily-a3pt5-evec3.gif", fps = 2)
    #    gif(anim, "ecmwf-15jul-20july-12hrly-a3pt5-evec3.gif", fps = 2)
    #gif(anim, "test.gif", fps=4)


end

function plot_spatemp_IDL(grid, V)

    N = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / N)

    K = size(V, 2)

    meanvariance = [mean([var(V[(t-1)*N+1:t*N, k]) for t = 1:T]) for k = 1:K]

    # Visualise the results

    x_disp = 1:K
    display(Plots.scatter(x_disp, meanvariance))

    return meanvariance

    #should be close to zero for temporal eigenfunctions (and the trivial eigenfunction) and far from zero for spatial eigenfunctions

end

function plot_gif(V, grid, vecnum)

    spacelength = length(grid.lonrange) * length(grid.latrange)
    T = Int(size(V)[1] / spacelength)
    xy = [[x, y] for t = 1:T for x ∈ grid.lonrange for y ∈ grid.latrange]
    xy = hcat(xy...)'

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    anim = @animate for t = 1:T
        display(heatmap(grid.lonrange, grid.latrange, reshape(real.(sliceV[t][:, vecnum]), length(grid.latrange), length(grid.lonrange)), c=:RdBu, xlabel="lon", ylabel="lat", title="First Spatial Eigenfunction (v_$vecnum), Day $(div(t-1,2)), $(((t-1)%2)*12):00"))
    end

    gif(anim, "InfDL_AlmostInvSets_a1p82_19Dec23.gif", fps=2)

end

#=
####annoying code if I want to find the lon/lat coordinates corresponding to an integer indexing G
Gindex = 1
dindex = findall(values(d) .== Gindex)
collect(d)[dindex]
=#
=#