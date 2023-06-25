using StatsBase
using LinearAlgebra

"""
    gaussian_wasserstein_metric(μ1, Σ1, μ2, Σ2)
    gaussian_wasserstein_metric(points1, points2)

The Wasserstein metric between two multivariate Gaussians `X_1 ~ N(μ1, Σ1)`
and `X_2 ~ N(μ2, Σ2)` is

    d^2 = ||μ1 - μ2||^2 + tr(Σ1 + Σ2 - 2*sqrt(Σ1*Σ2))

Some versions use `covAvg = sqrt(sqrt(Σ1) * Σ2 * sqrt(Σ1))`. 
This will not change `d^2`.
This is because `sqrt(Σ1) * Σ2 * sqrt(Σ1)` and `Σ1*Σ2` are similar matrices.
Hence they have the same eigenvalues and hence the same trace.

Source: https://en.wikipedia.org/wiki/Wasserstein_metric#Normal_distributions
"""
function gaussian_wasserstein_metric(μ1::AbstractMatrix, Σ1::AbstractMatrix, μ2::AbstractMatrix, Σ2::AbstractMatrix
    ; atol::AbstractFloat=1e-3)
    diff = μ1 - μ2
    covAvg = sqrt(Σ1 * Σ2)
    if eltype(covAvg) <: Complex
        @warn("sqrt(Σ1 * Σ2) is complex")
        if all(isapprox.(0.0, diag(imag(covAvg)), atol=atol))
            @info("imaginary components are small and have been set to zero")
            covAvg = real(covAvg)
        end
    end
    sum(diff .* diff) + tr(Σ1 + Σ2 - 2 * covAvg)
end

function gaussian_wasserstein_metric(points1::AbstractMatrix, points2::AbstractMatrix; options...)
    μ1 = mean(points1, dims=2)
    Σ1 = cov(points1, dims=2)
    μ2 = mean(points2, dims=2)
    Σ2 = cov(points2, dims=2)
    gaussian_wasserstein_metric(μ1, Σ1, μ2, Σ2; options...)
end
