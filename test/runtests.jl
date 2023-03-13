using RandomWalkBVP
using Test

#using Tullio
#using Folds
using Random
#using HDF5

mat = h5read("hdf5mat", "placa") |> BitArray;
@testset "Construction and basic functions" begin

    #@test valid(mat, 200, 200)

    #minimat = [false false false; false true false; false false false;] |> BitArray;
    minimat = [0 0 0; 0 1 0; 0 0 0];
    f(y, x) = 1 / (y + 1)
    RandomWalk(minimat, f) isa RandomWalk
    Random.seed!(42)
    sol = zeros(size(minimat))
    walkers = zeros(size(minimat))
    RandomWalk(minimat)
    @test sum(walkers) == 0
    @test sum(sol) == 0
    @test sum(minimat) == 1
    sol = trayectoria!(2,2, sol, minimat, walkers, f)
    @test sum(walkers) == 0
    @show sol
    @test sum(sol) ≈ 1/3
    @test sum(minimat) == 1
    for _ in 1:100
        @test walk!(1,1) ∈ ((2,1), (0, 1), (1, 2), (1, 0))
    end
    @test 0 == @allocated walk!(1,1)
end


@testet "Shrink grid" begin
    shs = shrink_grid(mat);
    @test size(shs) == (335, 389)
    @test 2 == findfirst(any.(>(0), eachrow(shs)))
    @test 334 == findlast(any.(>(0), eachrow(shs)))
    @test 2 == findfirst(any.(>(0), eachcol(shs)))
    @test 388 == findlast(any.(>(0), eachcol(shs)))
    # TODO idempotencia shs == eng
    @test shs == shrink_grid(shs)
end

# INPUT: 
# mat - BitMatrix de puntos interiores
function resuelve_Laplace2(mat, f, n)
    mat = shrink_grid(mat)
    n1, n2 = size(mat)
    sol = zeros(Float64, (n1, n2))
    walkers = zeros(Float64, (n1, n2))

    # TODO - reuse a single pushed_vector buffer
    # setup PushVector{(T,T)}() 
    view_walkers = @view walkers[mat]

    while all(<(n), view_walkers)
        for j in 1:n2
            for i in 1:n1
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
    return sol, numero_de_caminantes
end

@time sol, num = resuelve_Laplace2(mat, f, 1);