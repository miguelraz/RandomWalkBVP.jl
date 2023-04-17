using RandomWalkBVP
using HDF5
using Random
using Test

#using Tullio
#using Folds
#using Random

#@testset "RandomWalkBVP" showtiming=true begin
@testset "Construction" showtiming = true begin
    mat = h5read("hdf5mat", "placa")
    f(x, y) = 1 / (y + 1)
    re = RandomEnsemble(mat, f)
    @test all(==(0), re.temp_walkers)
    @test all(==(0), re.sol)
    @test all(==(0), re.walkers)
end
@testset "Basic functions" showtiming = true begin

    m = [0 0 0; 0 1 0; 0 0 0]
    b = BitArray(first(m) .!= m)
    f(x, y) = 1 / (y + 1)
    w = zeros(3, 3)
    tw = zeros(3, 3)
    sol = zeros(3, 3)
    vb = @view m[b]
    re = RandomEnsemble{Float64,typeof(f)}(f, b, w, tw, sol) #=vb, =#
    Random.seed!(42)
    @test re.bitmat == [false false false; false true false; false false false]
    @test sum(re.walkers) == 0
    # Only clear before using!
    @test sum(re.temp_walkers) == 0
    @test sum(re.sol) == 0
    @test sum(re.bitmat) == 1
    trajectory!(re, 2, 2)
    # 
    @test sum(re.temp_walkers) == 0
    @test sum(re.walkers) == 1
    @test sum(re.sol) ≈ 0.5
    @test sum(re.bitmat) == 1

    for _ in 1:100
        @test walk!(1, 1) ∈ ((2, 1), (0, 1), (1, 2), (1, 0))
    end
    # TODO - Fix ... ?
    #@test @allocations(walk!(1,1)) == 0
end


@testset "Shrink grid" showtiming = true begin
    mat = h5read("hdf5mat", "placa")
    shs = shrink_grid(mat)
    @test size(shs) == (335, 389)
    @test 2 == findfirst(any.(>(0), eachrow(shs)))
    @test 334 == findlast(any.(>(0), eachrow(shs)))
    @test 2 == findfirst(any.(>(0), eachcol(shs)))
    @test 388 == findlast(any.(>(0), eachcol(shs)))
    # TODO idempotencia shs == eng
    @test shs == shrink_grid(shs)
end

@testset "solve Laplace 3x3 single step" showtiming = true begin
    m = [0 0 0; 0 1 0; 0 0 0]
    b = BitArray(first(m) .!= m)
    f(x, y) = 1 / (y + 1)
    w = zeros(3, 3)
    tw = zeros(3, 3)
    sol = zeros(3, 3)
    vb = @view m[b]
    re = RandomEnsemble{Float64,typeof(f)}(f, b, w, tw, sol) #=vb,=#
    Random.seed!(42)

    @test sum(re.temp_walkers) == 0
    solve!(re, 1)
    @test sum(re.sol) ≈ 0.5
    @test sum(re.walkers) == 1
    @test sum(re.temp_walkers) == 0
end

@testset "Solve Laplace 334x389" showtiming = true begin
    using HDF5
    using Random
    mat = h5read("hdf5mat", "placa")
    f(x, y) = 1 / (y + 1)
    re = RandomEnsemble(mat, f)
    Random.seed!(42)
    solve!(re, 1)
end
#end