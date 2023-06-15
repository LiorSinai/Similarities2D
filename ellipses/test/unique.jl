using Test

@testset "no_repeat_columns" begin
    X = rand(10, 10)

    unique_idxs, counts = count_repeats(X)

    @test unique_idxs == 1:10
    @test counts == ones(10)
end;

@testset "repeat_columns_small" begin
    X = rand(2, 2)
    Y = X[:, [1, 2, 1]] + randn(2, 3) .* 1e-7

    unique_idxs, counts = count_repeats(Y)

    @test unique_idxs == [1, 2]
    @test counts == [2, 1]
end;

@testset "repeat_columns_large" begin
    n = 100
    X = rand(100, n)
    idxs = vcat([repeat([i], rand(1:5)) for i in 1:n]...)

    Y = X[:, idxs]
    Y += randn(size(Y)) .* 1e-7

    unique_idxs, counts = count_repeats(Y)
    Z = Y[:, unique_idxs]

    @test size(Z) == size(X)
    @test maximum(abs.(Z - X)) < 1e-6
end;