using StatsBase
using LinearAlgebra

"""
    gaussian_wasserstein_metric(μ1, Σ1, μ2, Σ2)
    gaussian_wasserstein_metric(points1, points2)

The Wasserstein metric between two multivariate Gaussians `X_1 ~ N(μ1, Σ1)`
and `X_2 ~ N(μ2, Σ2)` is
    d^2 = ||μ1 - μ2||^2 + tr(Σ1 + Σ2 - 2*sqrt(Σ1*Σ2))

Source: https://en.wikipedia.org/wiki/Wasserstein_metric#Normal_distributions
"""
function gaussian_wasserstein_metric(μ1::AbstractMatrix, Σ1::AbstractMatrix, μ2::AbstractMatrix, Σ2::AbstractMatrix)
    diff = μ1 - μ2
    cov2root = sqrt(Σ2) 
    covmean = sqrt(cov2root * Σ1 * cov2root)
    if eltype(covmean) <: Complex
        @warn("sqrt(Σ1 * Σ2) is complex")
        if all(isapprox.(0.0, diag(imag(covmean)), atol=1e-3))
            @info("imaginary components are small and have been set to zero")
            covmean = real(covmean)
        end
    end
    sum(diff .* diff) + tr(Σ1 + Σ2 - 2 * covmean)
end

function gaussian_wasserstein_metric(points1::AbstractMatrix, points2::AbstractMatrix)
    μ1 = mean(points1, dims=2)
    Σ1 = cov(points1, dims=2)
    μ2 = mean(points2, dims=2)
    Σ2 = cov(points2, dims=2)
    gaussian_wasserstein_metric(μ1, Σ1, μ2, Σ2)
end
