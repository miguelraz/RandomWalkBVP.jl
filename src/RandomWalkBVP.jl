module RandomWalkBVP
using Tullio
using LoopVectorization

__precompile__(false)
abstract type AbstractRandomWalkerBVP end
mutable struct RandomEnsemble{T,F} <: AbstractRandomWalkerBVP
    const f::F
    # TODO - make these subarrays
    const bitmat::BitMatrix
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

# TODO - Matrix{Bool} hack :(
function RandomEnsemble(xs::Matrix{T}, f::F) where {T,F}
    if T == Bool
        TT = Float64
    else
        TT = T
    end
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
    n, m = size(re.sol)
    if !valid(re, i, j)
        return
    end
    #@assert all(==(0), re.temp_walkers)::Bool
    myzero = zero(eltype(T))
    #re.temp_walkers .= zero(eltype(T))
    #@tullio re.temp_walkers[i, j] = myzero
    #@turbo for j in 1:m
    #    for i in 1:n
    #        re.temp_walkers[i, j] = myzero
    #    end
    #end

    # TODO -> keep top/bot/left/right maximums and extract that window
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

    # update!
    #re.sol[:,:] .+= re.temp_walkers .* scalar
    # Idea - single pass?
    #@tullio re.sol[i, j] += re.temp_walkers[i, j] * scalar
    #re.walkers .+= re.temp_walkers
    #@tullio re.walkers[i, j] += re.temp_walkers[i, j]
    @turbo for j in left:right
        for i in top:bot
            re.walkers[i, j] += re.temp_walkers[i, j]
            re.sol[i, j] += re.temp_walkers[i, j] * scalar
            re.temp_walkers[i, j] = myzero
        end
    end

    #@assert any(!=(0), re.temp_walkers)::Bool
    # Actually only need to clear before you use it, not twice!
    #re.temp_walkers .= zero(eltype(T))
    #@assert all(==(0), re.temp_walkers)::Bool
    return
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
export valid
export walk!
export trajectory!
export shrink_grid
export solve!
export RandomWalkBVP

end
