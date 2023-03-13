module RandomWalkBVP
__precompile__(false)
abstract type AbstractRandomWalkerBVP end
mutable struct RandomEnsemble{T} <: AbstractRandomWalkerBVP
    const f::Function
    # TODO - make these subarrays
    const bitmat::AbstractArray
    const view_bitmat::AbstractArray
    walkers::Matrix{T}
    temp_walkers::Matrix{T}
    sol::Matrix{T}
end

function RandomEnsemble(xs::Matrix{T}, f) where {T}
    n, m = size(xs)
    mat = first(xs) .== xs
    b = shrink_grid(mat)
    vb = @view b[b]
    w = zeros(T, n, m)
    tw = zeros(T, n, m)
    s = zeros(T, n, m)
    @show typeof(f)
    @show typeof(b)
    @show typeof(vb)
    @show typeof(w)
    @show typeof(tw)
    @show typeof(s)
    RandomEnsemble{T}(f, b, vb, w, tw, s)
end

# TODO 1-liner?
function walk!(re::RandomEnsemble{T}, i, j) where {T}
    n = rand(((0, 1), (0, -1), (1, 0), (-1, 0)))
    return i + n[1], j + n[2]
end
valid(re::RandomEnsemble{T}, i, j) where {T} = re.bitmat[i, j]

function trajectory!(re::RandomEnsemble{T}, i, j) where {T}
    @assert all(==(zero(eltype(T))), re.temp_walkers)
    re.temp_walkers .= zero(eltype(T))

    # TODO -> keep top/bot/left/right maximums and extract that window
    while valid(re, i, j)
        i, j = walk!(re, i, j)
        re.temp_walkers[i, j] += 1
    end
    scalar = re.f(i, j)

    re.sol .+= re.temp_walkers .* scalar

    @assert re.temp_walkers != zeros(size(re.temp_walkers))
    re.temp_walkers .= zero(eltype(T))
    @assert all(==(0), re.temp_walkers)
    return re.sol
end

function shrink_grid(grid)
    right = findlast(any.(>(0), eachcol(grid)))
    bot = findlast(any.(>(0), eachrow(grid)))
    any(<(3), (right, bot)) && error("Grid too small")
    # TODO - add warnings for short boundaries on right/bot
    n, m = size(grid)
    left = findfirst(any.(>(0), eachcol(grid)))
    top = findfirst(any.(>(0), eachrow(grid)))

    # 
    left -= left != 1
    top -= top != 1
    bot += bot != m
    right += right != n

    @view grid[top:bot, left:right]
end

function solve(re::RandomEnsemble{T}, n::Int) where {T}
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

export RandomEnsemble
export valid
export walk!
export trajectory!
export shrink_grid

end
