using Random
using StatsBase
using LinearAlgebra
using Printf
using PolygonAlgorithms: matrix_to_points

include("wasserstein.jl")
include("similarity_metrics.jl")
include("bivariate.jl")

function generate_distributions(rng::AbstractRNG, r0::T, μ0::T, σ0::T, sizes::UnitRange{Int};
    N::Int=100) where T <: AbstractFloat
    r = r0 * rand(rng, N) .- 1
    μ = μ0 * (rand(rng, 2, N) .- 0.5)
    σ = σ0 * rand(rng, 2, N) .+ 0.9
    n = rand(rng, sizes, N)
    r, μ, σ, n
end

num_distributions = 1000
rng = MersenneTwister(15)

rs, μs, σs, ns = generate_distributions(rng, 2.0, 10.0, 6.1, 1_000:10_000; N=num_distributions);

idx = 1
x1, y1 = correlated_bivariate_normal(μs[:, idx], σs[:, idx], rs[idx], ns[idx])
points1 = permutedims(hcat(x1, y1))
points1_vec = matrix_to_points(points1)

stats = Dict(
    :wasserstein => zeros(num_distributions),
    :ellipse => zeros(num_distributions),
    :polygon => zeros(num_distributions),
    :polygon2 => zeros(num_distributions),
)
for j in 1:num_distributions
    x2, y2 = correlated_bivariate_normal(μs[:, j], σs[:, j], rs[j], ns[j])
    points2 = permutedims(hcat(x2, y2))

    points2_vec = matrix_to_points(points2)
    # initialise
    similarities_ellipses(points1, points2)
    similarities_polygons(points1, points2)
    gaussian_wasserstein_metric(points1, points2)

    stats_ellipse = @timed similarities_ellipses(points1, points2);
    stats_polygon = @timed similarities_polygons(points1_vec, points2_vec);
    stats_wasserstein = @timed gaussian_wasserstein_metric(points1, points2);
    stats[:polygon][j] += stats_polygon.time
    stats[:ellipse][j] += stats_ellipse.time
    stats[:wasserstein][j] += stats_wasserstein.time
end

base = stats[:wasserstein]
time_ellipses = stats[:ellipse] ./ base
time_polygons = stats[:polygon] ./ base

@printf("Times relative to Wasserstein normal approximation:\n")
@printf("%10s      min       max      mean ± std   \n", "")
for (method, times) in zip(["ellipses", "polygons"], [time_ellipses, time_polygons])
    @printf("%10s %8.4f  %8.4f  %8.4f ± %.4f\n", method, minimum(times), maximum(times), mean(times), std(times))
end
