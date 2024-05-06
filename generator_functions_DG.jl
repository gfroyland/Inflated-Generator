using Plots, LinearAlgebra, QuadGK, SparseArrays, Arpack, Statistics, ProgressMeter

#create a data structure for the grid;  it might be expanded later to include e.g. G, spectrum, evecs, etc... or other items
struct Grid
    centres
    x_range
    y_range
    x_min
    x_max
    y_min
    y_max
    Δ_x
    Δ_y
end

include("SEBA.jl")

#create a dictionary to do the indexing we want and a grid struct
function make_dict_grid(x_min, x_max, Δ_x, y_min, y_max, Δ_y)
    x_range = x_min+Δ_x/2:Δ_x:x_max-Δ_x/2
    y_range = y_min+Δ_y/2:Δ_y:y_max-Δ_y/2
    x_range = round.(x_range, digits=6)
    y_range = round.(y_range, digits=6)
    i = 0
    temparray = []
    for x ∈ x_range
        for y ∈ y_range
            i += 1
            push!(temparray, ([x, y], i))
        end
    end
    d = Dict(temparray)
    centres = Tuple.(collect(keys(d)))
    grid = Grid(centres, x_range, y_range, x_min, x_max, y_min, y_max, Δ_x, Δ_y)
    return d, grid
end

function make_generator(d, grid, F, ϵ)

    #create list of box centres
    centres = collect(keys(d))
    #compute volumes of cells. constant at the moment, but adjust later
    #WHEN MAKING IT NONCONSTANT LATER, BE SURE TO MATCH THE CORRECT INDEXING IN THE DICTIONARY AND ENSURE IT IS INDEXED CORRECTLY IN THE G COMPUTATION
    #GF: I believe it is fixed below.

    volume = zeros(length(centres))
    for c ∈ centres
        volume[d[c]] = grid.Δ_x * grid.Δ_y
    end

    #create basic lon,lat increment vectors to access adjacent grid cells
    δx = [grid.Δ_x, 0]
    δy = [0, grid.Δ_y]

    #create list of box centres
    #create an array G (for Generator) to hold the flux values. 
    #G is the main array we want to compute
    x_length = length(grid.x_range)
    y_length = length(grid.y_range)
    Gdim = x_length * y_length
    G = spzeros(Gdim, Gdim)

    #set up normal vectors to each of the 4 faces and a parametric function for each face
    rightnormal = δx / norm(δx)  #normal vector pointing to the right...similarly make another 3, one for each direction
    rightface(c, t) = c + (δx / 2 + δy / 2) - δy * t   #parameterise right face of box with centre t by a parameter t ranging from 0 to 1
    leftnormal = -δx / norm(δx)
    leftface(c, t) = c + (-δx / 2 + δy / 2) - δy * t
    uppernormal = δy / norm(δy)
    upperface(c, t) = c + (δx / 2 + δy / 2) - δx * t
    lowernormal = -δy / norm(δy)
    lowerface(c, t) = c + (δx / 2 - δy / 2) - δx * t

    #construct the generator matrix G
    tol = 1e-2  #hard coded for now, but ultimately should be adapted to entry sizes of 
    intorder = 1
    #start looping over c (centres of cells)
    for c ∈ centres
        rightc = round.(c + δx, digits=6)
        if rightc ∈ keys(d)  #check that the box on the right exists
            #compute the entry for G corresponding to flux through the right face
            #I add ϵ directly to the flux value, so it is also integrated from 0 to 1
            #Because of the integration from 0 to 1 instead of 0 to the length of the face, instead of
            #dividing by volume, I multiply by (Δ_y/volume), which is 1/Δ_x.
            #similarly in the other faces below
            #G[d[c], d[rightc]] = (quadgk(t -> max(F(rightface(c, t)) ⋅ rightnormal, 0) + ϵ, 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_x
            G[d[c], d[rightc]] = (quadgk(t -> max(F(rightface(c, t)) ⋅ rightnormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_x + (ϵ^2 / (2 * (grid.Δ_x)^2))
        end
        leftc = round.(c - δx, digits=6)
        if leftc ∈ keys(d)
            #G[d[c], d[leftc]] = (quadgk(t -> max(F(leftface(c, t)) ⋅ leftnormal, 0) + ϵ, 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_x
            G[d[c], d[leftc]] = (quadgk(t -> max(F(leftface(c, t)) ⋅ leftnormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_x + (ϵ^2 / (2 * (grid.Δ_x)^2))
        end
        upperc = round.(c + δy, digits=6)
        if upperc ∈ keys(d)
            #G[d[c], d[upperc]] = (quadgk(t -> max(F(upperface(c, t)) ⋅ uppernormal, 0) + ϵ, 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_y
            G[d[c], d[upperc]] = (quadgk(t -> max(F(upperface(c, t)) ⋅ uppernormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_y + (ϵ^2 / (2 * (grid.Δ_y)^2))
        end
        lowerc = round.(c - δy, digits=6)
        if lowerc ∈ keys(d)
            #G[d[c], d[lowerc]] = (quadgk(t -> max(F(lowerface(c, t)) ⋅ lowernormal, 0) + ϵ, 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_y
            G[d[c], d[lowerc]] = (quadgk(t -> max(F(lowerface(c, t)) ⋅ lowernormal, 0), 0, 1, rtol=tol, atol=tol, order=intorder)[1]) / grid.Δ_y + (ϵ^2 / (2 * (grid.Δ_y)^2))
        end
    end
    #G = G - spdiagm(vec(sum(G, dims=2)))

    #place negative row sums on the diagonal of G so that the row sum of G is now zero.
    
    G = G - spdiagm(vec(sum(spdiagm(1 ./ volume) * G * spdiagm(volume), dims=2)))

    #adjust G by a similarity transformation to ensure that the matrix has row sum 0
    #this ensures leading right evec is constant and leading left evec is a *density* rather than a measure on each box.
    #see Lemma 4.7 in FJK'13.

    G = spdiagm(1 ./ volume) * G * spdiagm(volume)

    #I suspect I can replace this lengthy code with the code above...indexing same?
    #=
    Gdiag = spzeros(Gdim, Gdim)
    for c1 ∈ centres
        for c2 ∈ centres
            Gdiag[d[c1], d[c2]] = volume[d[c2]] / volume[d[c1]] * G[d[c1], d[c2]]
        end
    end
    G = G - spdiagm(vec(sum(Gdiag, dims=2)))
    =#

    return G
end

function make_dynamic_generator(Gvec)

    Gᴰ = mean(Gvec)
    return Gᴰ

end

function make_inflated_generator(Gvec, Δt, a)

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
    𝐋 = kron(L/(Δt)^2, one(Gvec[1]))     #second argument is identity matrix of size Gvec[1]
    #additively combine Gspat and Gtemp to create inflated generator
    𝐆 = 𝐆spat + a^2 * 𝐋

    return 𝐆
end

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
    FF(x, y) = F([x, y])
    #display(surface(grid.lonrange, grid.latrange, norm ∘ FF, xlabel="lon", ylabel="lat", zlabel="wind speed (units?)", c=:jet))
    display(Plots.contourf(grid.x_range, grid.y_range, norm ∘ FF, xlabel="x", ylabel="y", title="Double Gyre Speed", aspectratio=1, c=:jet))

    #plot quiver of windspeed field, with scaling of arrows inspired by the matlab autoscaling
    centres = collect(keys(d))
    sparsity = 7  #include only every sparsity^{th} grid point, ideally sparsity is coprime with the number of grid points in either lon or lat direction
    sparsecentres = sort(centres)[1:sparsity:end]
    sparsecentremat = hcat(sparsecentres...)'
    #vf=F.(sparsecentres)
    #vfmat=vcat(vf...)
    meannorm = mean(norm.(F.(sparsecentres)))  #mean vector norm
    scale = norm([grid.Δ_x, grid.Δ_y]) * (√sparsity) / meannorm #scale quiver vectors by multiplicative factor "scale", computed so that the average-length vector doesn't overlap
    vF(x, y) = scale * FF(x, y)
    display(quiver(sparsecentremat[:, 1], sparsecentremat[:, 2], quiver=vF, title="DG Vector Field", arrow=:filled, aspectratio=1, xlabel="x", ylabel="y"))

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
    p = []
    for k = 1:9
        push!(p, heatmap(grid.x_range, grid.y_range, real.(reshape(v[:, k], length(grid.y_range), length(grid.x_range))), c=:RdBu, title="eigenvector $k", aspectratio=1))
    end
    display(Plots.plot(p..., layout=(3, 3)))

    #plots leading seba vectors
    pseba = []
    for k = 1:size(s)[2]
        push!(pseba, heatmap(grid.x_range, grid.y_range, real.(reshape(s[:, k], length(grid.y_range), length(grid.x_range))), c=:Reds, title="SEBA vector $k", aspectratio=1))
    end
    display(Plots.plot(pseba..., layout=(3, 3)))

    #plot max of seba vectors
    smax = maximum(s, dims=2)
    #display(heatmap(grid.lonrange, grid.latrange, real.(reshape(smax, length(grid.latrange), length(grid.lonrange))), title="superposition of $(size(s)[2]) SEBA vectors", c=:Reds, aspectratio=1, xlabel="lon", ylabel="lat"))

    #display(real.(reshape(smax, length(grid.latrange), length(grid.lonrange))), grid.lonmin:grid.lonspacing:grid.lonmax,  grid.latmin:grid.latspacing:grid.latmax, projection=:Robinson, coast=true, colorbar=true, cmap=:red2green, fmt=:png, title="superposition of $(size(s)[2]) SEBA vectors")
    heatmap(grid.x_range, grid.y_range, real.(reshape(smax, length(grid.y_range), length(grid.x_range))), title="superposition of $(size(s)[2]) SEBA vectors", c=:jet, aspectratio=1, xlabel="x", ylabel="y")

end


function plot_slices(V, grid, vecnum)

    #X=[(t,x,y) for t=1:T  for x∈grid.lonrange for y∈grid.latrange ]
    #scatter(X,zcolor=real(V[:,2]))

    spacelength = length(grid.x_range) * length(grid.y_range)
    T = Int(size(V)[1] / spacelength)
    xy = [[x, y] for t = 1:T for x ∈ grid.x_range for y ∈ grid.y_range]
    xy = hcat(xy...)'

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]
    slicenorm = [norm(sliceV[t][:, k]) for t = 1:T, k = 1:min(5, size(V)[2])]
    slicemean = [mean(sliceV[t][:, k]) for t = 1:T, k = 1:min(5, size(V)[2])]
    display(Plots.plot(Plots.plot(real(slicenorm)), Plots.plot(real(slicemean)), layout=(2, 1)))

    cmax = maximum(abs.(real(V[:, vecnum])))
    anim = @animate for t = 1:T
        tm = (t - 1) / (T - 1)
        display(heatmap(grid.x_range, grid.y_range, reshape(real.(sliceV[t][:, vecnum]), length(grid.y_range), length(grid.x_range)), clims=(-cmax,cmax), c=:RdBu, xlabel="x", ylabel="y", title="t = $tm"))
        #display(imshow(reshape(real.(sliceV[t][:, vecnum]), length(grid.latrange), length(grid.lonrange)), grid.lonmin:grid.lonspacing:grid.lonmax,  grid.latmin:grid.latspacing:grid.latmax, projection=:Robinson, coast=true, colorbar=true, cmap=:red2green, fmt=:png, title="Day $t"))

        #imshow(-G, -5:0.5:42, 45:0.5:75, projection=:Robinson, coast=true, colorbar=true, fmt=:png, cmap=:red2green)
        #also cmap: panoply, hot, bam, vik, balance, buda, haline, matter, amp, oslo, gebco, curl
    end

    #   gif(anim, "ecmwf-24jul-4aug-daily-a3pt5-evec3.gif", fps = 2)
    #    gif(anim, "ecmwf-15jul-20july-12hrly-a3pt5-evec3.gif", fps = 2)
    gif(anim, "DG_inflated_slices_vector$vecnum.mp4", fps=4)


end

function plot_9vecs_IDL(grid, Λ, V)

    #plot inflated spectrum
    display(Plots.scatter(Λ, title="$(length(Λ)) eigenvalues with largest real part, a = $a"))

    spacelength = length(grid.x_range) * length(grid.y_range)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    cmax = maximum(abs.(real(V)), dims=1)

    for τ ∈ 1:T
    #anim = @animate for τ ∈ 1:T

        τ_m = (τ - 1) / (T - 1)
        P = []
        for k = 1:9
            push!(P, Plots.contourf(grid.x_range, grid.y_range, reshape(real.(sliceV[τ][:, k]), length(grid.y_range), length(grid.x_range)), clims=(-cmax[k], cmax[k]), c=:RdBu, xlabel="x", ylabel="y", title="v_$k, t = $τ_m", linewidth=0, levels=100))
        end
        display(Plots.plot(P..., layout=(3, 3)))
    end

    #gif(anim, "EuroBlock_0DayExt_9EigVecs_a3p5_30Jan24.gif", fps=2)
    #gif(anim, "DoubleGyre_9EigVecs_22Apr24_2.mp4", fps=2)

end

function plot_spatemp_IDL(grid, Λ, V)

    #plot inflated spectrum
    display(Plots.scatter(Λ, title="$(length(Λ)) eigenvalues with largest real part, a = $a"))

    N = length(grid.x_range) * length(grid.y_range)
    T = Int(size(V)[1] / N)

    K = size(V, 2)

    meanvariance = [mean([var(V[(t-1)*N+1:t*N, k]) for t = 1:T]) for k = 1:K]

    # Visualise the results

    x_disp = 1:K
    display(Plots.scatter(x_disp, meanvariance))

    return meanvariance

    #should be close to zero for temporal eigenfunctions (and the trivial eigenfunction) and far from zero for spatial eigenfunctions

end

function plot_SEBA_IDL(grid, Σ)

    # Prepare the slices of nine eigenvectors at each time step 
    spacelength = length(grid.x_range) * length(grid.y_range)
    T = Int(size(Σ)[1] / spacelength)

    sliceΣ = [Σ[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    # Plot All SEBA Evolutions

    for σ ∈ 1:size(Σ)[2]
        vidname = "DG_SEBA_3Vecs_" * lpad(σ, 1, "0") * "_a1p5_31Jan24.mp4"
        for t ∈ 1:T
        #anim = @animate for t ∈ 1:T

            τ_m = (t - 1) / (T - 1)
            display(Plots.contourf(grid.x_range, grid.y_range, reshape(sliceΣ[t][:, σ], length(grid.y_range), length(grid.x_range)), clims=(0, 1), c=:Reds, xlabel="x", ylabel="y", title="SEBA $σ, t = $τ_m", linewidth=0, levels=100))

        end
        #gif(anim, vidname, fps=2)
        println("SEBA Video $σ Complete.")
    end

    # Plot the Maxima of each SEBA Vector

    vidname = "DG_SEBAMax_3Vecs_a1p5_31Jan24.mp4"
    for 𝒯 ∈ 1:T
    #anim = @animate for 𝒯 ∈ 1:T

        τ_m = (𝒯 - 1) / (T - 1)
        display(Plots.contourf(grid.x_range, grid.y_range, reshape(maximum(sliceΣ[𝒯][:, :], dims=2), length(grid.y_range), length(grid.x_range)), clims=(0, 1), c=:Reds, xlabel="x", ylabel="y", title="Max SEBA, t = $τ_m", linewidth=0, levels=100))

    end
    #gif(anim, vidname, fps=2)
    println("Max SEBA Video Complete.")

end

function save_results(grid, λ, v, s, Λ, V, Σ, filename)

    file_ID = h5open(filename, "w")

    file_ID["longitude"] = grid.x_range
    file_ID["latitude"] = grid.y_range

    s_max = maximum(s, dims=2)
    Σ_max = maximum(Σ, dims=2)

    file_ID["DL_Lambda_Real"] = real.(λ)
    file_ID["DL_Lambda_Imag"] = imag.(λ)
    for β ∈ 1:size(v, 2)
        varname = "DL_v_" * lpad(β, 2, "0")
        file_ID[varname] = real.(reshape(v[:, β], length(grid.y_range), length(grid.x_range)))
    end
    for γ ∈ 1:size(s,2)
        varname = "DL_SEBA_" * lpad(γ, 2, "0")
        file_ID[varname] = real.(reshape(s[:, γ], length(grid.y_range), length(grid.x_range)))
    end
    file_ID["DL_SEBA_Max"] = real.(reshape(s_max, length(grid.y_range), length(grid.x_range)))

    file_ID["IDL_Lambda_Real"] = real.(Λ)
    file_ID["IDL_Lambda_Imag"] = imag.(Λ)

    spacelength = length(grid.x_range) * length(grid.y_range)
    T = Int(size(V)[1] / spacelength)

    sliceV = [V[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    sliceΣ = [Σ[(t-1)*spacelength.+(1:spacelength), :] for t = 1:T]

    for τ ∈ 1:T

        for φ ∈ 1:size(V, 2)
            varname = "IDL_v_" * lpad(φ, 2, "0") * "_t" * lpad(τ, 2, "0")
            file_ID[varname] = real.(reshape(sliceV[τ][:, φ], length(grid.y_range), length(grid.x_range)))
        end
        for ψ ∈ 1:size(Σ, 2)
            varname = "IDL_SEBA_" * lpad(ψ, 2, "0") * "_t" * lpad(τ, 2, "0")
            file_ID[varname] = real.(reshape(sliceΣ[τ][:, ψ], length(grid.y_range), length(grid.x_range)))
        end

        varname = "IDL_SEBA_Max_t" * lpad(τ, 2, "0")
        file_ID[varname] = real.(reshape(Σ_max[(τ-1)*spacelength.+(1:spacelength)], length(grid.y_range), length(grid.x_range)))
    end

    close(file_ID)

end

#=
####annoying code if I want to find the lon/lat coordinates corresponding to an integer indexing G
Gindex = 1
dindex = findall(values(d) .== Gindex)
collect(d)[dindex]
=#
