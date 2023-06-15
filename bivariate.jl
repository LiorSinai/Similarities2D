#https://www.probabilitycourse.com/chapter5/5_3_2_bivariate_normal_dist.php

function correlated_bivariate_normal(μ::Tuple{T, T}, σ::Tuple{T, T}, r::T, n::Int) where T <: AbstractFloat
    z1 = randn(n);
    z2 = randn(n);
    x = μ[1]  .+ z1 .* σ[1];
    y = μ[2] .+ σ[2] .* (r .* z1 .+ sqrt(1 - r^2) .* z2);
    return (x, y)
end

function correlated_bivariate_normal(μ::Vector{T}, σ::Vector{T}, r::T, n::Int) where T <: AbstractFloat
    z1 = randn(n);
    z2 = randn(n);
    x = μ[1]  .+ z1 .* σ[1];
    y = μ[2] .+ σ[2] .* (r .* z1 .+ sqrt(1 - r^2) .* z2);
    return (x, y)
end
