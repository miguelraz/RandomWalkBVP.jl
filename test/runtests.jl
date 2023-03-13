using RandomWalkBVP
using Test

#using Tullio
#using Folds
using Random
using HDF5

mat = h5read("hdf5mat", "placa") |> BitArray;
@testset "Construction and basic functions" begin
    #@test valid(mat, 200, 200)

    m = [0 0 0; 0 1 0; 0 0 0];
    b = BitArray(first(m) .!=  m)
    f(x, y) = 1 / (y + 1)
    w = zeros(3,3)
    tw = zeros(3,3)
    sol = zeros(3,3)
    vb = @view m[b]
    re = RandomEnsemble{Float64}(f, b, vb, w, tw, sol)
    Random.seed!(42)
    @test re.bitmat == [false false false; false true false; false false false]
    @test sum(re.walkers) == 0
    @test sum(re.sol) == 0
    @test sum(re.bitmat) == 1
    trajectory!(re, 2, 2)
    @test sum(re.walkers) == 0
    #@test sum(re.sol) ≈ 1 / 3
    @test sum(re.bitmat) == 1

    for _ in 1:100
        re = RandomEnsemble{Float64}(f, b, vb, w, tw, sol)
        @test walk!(re, 1, 1) ∈ ((2, 1), (0, 1), (1, 2), (1, 0))
    end
    # TODO - Fix ... ?
    #@test 0 == @allocated walk!(re, 1, 1)
end


@testset "Shrink grid" begin
    shs = shrink_grid(mat)
    @test size(shs) == (335, 389)
    @test 2 == findfirst(any.(>(0), eachrow(shs)))
    @test 334 == findlast(any.(>(0), eachrow(shs)))
    @test 2 == findfirst(any.(>(0), eachcol(shs)))
    @test 388 == findlast(any.(>(0), eachcol(shs)))
    # TODO idempotencia shs == eng
    @test shs == shrink_grid(shs)
end

@testset "solve Laplace" begin
    @test true
end
