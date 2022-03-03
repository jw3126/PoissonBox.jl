using PoissonBox
using Test
using IterativeSolvers
import PoissonBox as PB
using Random: MersenneTwister


@testset "laplace eqn operator" begin
    forward = PoissonBox.forward
    @testset "1d" begin
        @test forward(1:5) ≈ [1,0,0,0,5]
        @test forward([0,0,1,0,0]) ≈ [0,0.5,-1,0.5,0]
    end
    @testset "2d" begin
        @test forward(
            [0 0 0 0 0
                 0 0 0 0 0
                 0 0 1 0 0
                 0 0 0 0 0 
                 0 0 0 0 0 
                ]) ≈ [
                 0 0    0    0    0
                 0 0    0.25 0    0
                 0 0.25 -1   0.25 0
                 0 0    0.25 0    0 
                 0 0    0    0    0 
                ]
    end
    @testset "3d" begin
        arr = fill(0, (5,5,5))
        arr[3,3,3] = 1
        expected = fill(0.0,(5,5,5))
        expected[3,3,3] = -1
        expected[2,3,3] = 1/6
        expected[4,3,3] = 1/6
        expected[3,2,3] = 1/6
        expected[3,4,3] = 1/6
        expected[3,3,2] = 1/6
        expected[3,3,4] = 1/6
        @test forward(arr) ≈ expected
    end
    @testset "3d bdry" begin
        arr = randn(3,4,5)
        out = forward(arr)
        expected = copy(arr)
        expected[2:2, 2:3, 2:4] .= 0
        out[2:2, 2:3, 2:4] .= 0
        @test expected ≈ out
    end
end

@testset "solve" begin
    solve = PoissonBox.solve
    y = Float64[1, 0, 0, 0, 5]
    @test solve(y) ≈ 1:5
    for shape in [(5,),
                  (4,3),
                  (5,4),
                  (4,3,4),
                 ]
        y = rand(MersenneTwister(0), shape...)
        x = solve(y)
        @test PB.forward(x) ≈ y rtol=1e-3
    end
end

@testset "solve low level" begin
    op, v = PB.poisson_problem([1, 0, 0, 0, 5])
    sol = bicgstabl(op, v, 2)
    u = PB.reshape_solution(op, sol)
    @test u ≈ 1:5
end
