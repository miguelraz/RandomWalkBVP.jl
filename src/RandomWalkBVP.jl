module RandomWalkBVP

__precompile__(false)
abstract type AbstractRandomWalkerBVP end
mutable struct RandomEnsemble{T} <: AbstractRandomWalkerBVP
    const f::Function
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
    any(<(3), (right, bot))::Bool && error("Grid too small")

    # TODO - add warnings for short boundaries on right/bot
    n, m = size(grid)
    if n == right || m == bot
        error("Grid too small, add a right/bot edge")
    end

    left = findfirst(any.(>(0), eachcol(grid)))
    top = findfirst(any.(>(0), eachrow(grid)))

    @view grid[top-1:bot+1, left-1:right+1]
end

# TODO - Matrix{Bool} hack :(
function RandomEnsemble(xs::Matrix{T}, f) where {T}
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
    RandomEnsemble{TT}(f, b, #=vb,=# w, tw, s)
end

# TODO 1-liner?
@inline function walk!(re::RandomEnsemble{T}, i, j) where {T}
    n = rand(((0, 1), (0, -1), (1, 0), (-1, 0)))
    return i + n[1], j + n[2]
end
valid(re::RandomEnsemble{T}, i, j) where {T} = re.bitmat[i, j]

function trajectory!(re::RandomEnsemble{T}, i, j) where {T}
    if !valid(re, i, j)
        return re.sol
    end
    #@assert all(==(0), re.temp_walkers)::Bool
    re.temp_walkers .= zero(eltype(T))

    # TODO -> keep top/bot/left/right maximums and extract that window
    while valid(re, i, j)
        i, j = walk!(re, i, j)
        re.temp_walkers[i, j] += 1
    end
    scalar::T = re.f(i, j)

    # update!
    re.sol .+= re.temp_walkers .* scalar
    re.walkers .+= re.temp_walkers

    #@assert any(!=(0), re.temp_walkers)::Bool
    # Actually only need to clear before you use it, not twice!
    re.temp_walkers .= zero(eltype(T))
    #@assert all(==(0), re.temp_walkers)::Bool
    return re.sol
end


function solve!(re::RandomEnsemble{T}, k::Int) where {T}
    n, m = size(re.bitmat)
    # TODO - only iterate over valid indexes
    while all(<(k), re.walkers)::Bool
        for j in 1:m
            for i in 1:n
                trajectory!(re, i, j)
            end
        end
        all(>=(k), re.walkers)::Bool && break
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
