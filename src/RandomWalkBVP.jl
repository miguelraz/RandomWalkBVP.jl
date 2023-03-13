module RandomWalkBVP

abstract type AbstractRandomWalkerBVP end
mutable struct RandomEnsemble{T} <: AbstractRandomWalkerBVP
    const f::Function
    const bitmat::BitMatrix
    const view_bitmat::BitMatrix
    walkers::Matrix{T}
    temp_walkers::Matrix{T}
    sol::Matrix{T}
end

function RandomEnsemble{T}(xs::Matrix{T}, f::Function) where T
    n, m = size(xs)
    b = [xs[i, j] == first(xs) for i in 1:n, j in 1:m] |> BitArray
    b = shrink_grid(b)
    vb = @view b[b]
    w = zeros(T, n, m)
    tw = zeros(T, n, m)
    s = zeros(T, n, m)
    RandomWalk{T}(f, b, vb, w, tw, s)
end

walk!(re::RandomEnsemble{T}, i, j) where T = (i, j) .+ rand(((0, 1), (0, -1), (1, 0), (-1, 0)))
valid(re::RandomEnsemble{T}, i, j) where T = re.bitmat[i, j]

function trajectory!(re::RandomEnsemble{T} i, j) where T
    @assert all(==(zero(eltype(T))), walkers)
    re.tempwalkers .= zero(eltype(T))

    # TODO -> keep top/bot/left/right maximums and extract that window
    while valid(re, i, j)
        i, j = walk!(re, i, j)
        re.temp_walkers[i,j] += 1
    end
    scalar = f(i, j)

    re.sol .+= re.temp_walkers .* scalar

    @assert walkers != zeros(size(sol))
    walkers .= zero(eltype(T))
    @assert all(==(0), walkers)
    return sol
end

function shrink_grid(grid::BitMatrix)
    left = findfirst(any.(>(0), eachcol(grid)))
    right = findlast(any.(>(0), eachcol(grid)))

    top = findfirst(any.(>(0), eachrow(grid)))
    bot = findlast(any.(>(0), eachrow(grid)))
    # TODO - handle corner cases where ±1 might not be inside the array
    @view grid[top-1:bot+1, left-1:right+1]
end

function solve(re::RandomEnsemble{T}, n::Int) where T
    n, m = size(re.bitmat)
    while all(<(n), view_bitmat)
        for j in 1:m
            for i in 1:n
                # Walkers
                # Evolucionar caminante hasta frontera 
                    # agregar a matriz de walkers
                    # guardar valor en frontera
                    # mul! sol <- fend * num_walks
                caminante2!(i, j, mat, sol, numero_de_caminantes, f)
            end
            all(>=(n), view_walkers) && break
        end
    end
    return re
end

export RandomWalk
export valid
export walk!
export trajectory!
export shrink_grid

end
