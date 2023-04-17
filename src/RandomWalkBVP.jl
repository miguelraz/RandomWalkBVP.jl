module RandomWalkBVP
using LoopVectorization

abstract type AbstractRandomWalkerBVP end
mutable struct RandomEnsemble{T,F} <: AbstractRandomWalkerBVP
    const f::F
    const bitmat::BitMatrix
    # TODO - make these subarrays
    #const view_bitmat::AbstractArray{S, 2}
    walkers::Matrix{T}
    temp_walkers::Matrix{T}
    sol::Matrix{T}
end

function shrink_grid(grid)
    right = findlast(any.(>(0), eachcol(grid)))
    bot = findlast(any.(>(0), eachrow(grid)))
    (right < 3) || (bot < 3) && error("Grid too small")

    # TODO - add warnings for short boundaries on right/bot
    n, m = size(grid)
    (n == right) || (m == bot) && error("Grid too small, add a right/bot edge")

    left = findfirst(any.(>(0), eachcol(grid)))
    top = findfirst(any.(>(0), eachrow(grid)))

    @view grid[top-1:bot+1, left-1:right+1]
end

function RandomEnsemble(xs::Matrix{T}, f::F) where {T,F}
    TT = T == Bool ? Float64 : T # Handle Matrix{Bool}
    mat = first(xs) .!= xs
    b = shrink_grid(mat) |> BitArray
    n, m = size(b)
    vb = @view b[b]
    w = zeros(TT, n, m)
    tw = zeros(TT, n, m)
    s = zeros(TT, n, m)
    RandomEnsemble{TT,F}(f, b, w, tw, s) #=vb,=#
end

# TODO 1-liner?
@inline function walk!(i, j)
    n = rand(((0, 1), (0, -1), (1, 0), (-1, 0)))
    @inbounds i += n[1]
    @inbounds j += n[2]
    return i, j
end

valid(re::RandomEnsemble{T,F}, i, j) where {T,F} = re.bitmat[i, j]

function trajectory!(re::RandomEnsemble{T,F}, i, j) where {T,F}
    !valid(re, i, j) && return
    n, m = size(re.sol)

    top = i
    bot = i
    left = j
    right = j
    while valid(re, i, j)
        i, j = walk!(i, j)
        top = min(i, top)
        bot = max(i, bot)
        right = max(j, right)
        left = min(j, left)
        @inbounds re.temp_walkers[i, j] += 1
    end

    scalar = re.f(i, j)
    myzero = zero(eltype(T))
    @turbo for j in left:right
        for i in top:bot
            re.walkers[i, j] += re.temp_walkers[i, j]
            re.sol[i, j] += re.temp_walkers[i, j] * scalar
            re.temp_walkers[i, j] = myzero
        end
    end
end


function solve!(re::RandomEnsemble{T,F}, k::Int) where {T,F}
    n, m = size(re.bitmat)
    # TODO - only iterate over valid indexes
    while all(<(k), re.walkers)::Bool
        for j in 1:m
            for i in 1:n
                trajectory!(re, i, j)
            end
            all(>=(k), re.walkers)::Bool && break
        end
    end
    return re
end

export RandomEnsemble
export valid, walk!, trajectory!, solve!, shrink_grid
end