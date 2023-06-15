using Test

# basic expansion of the polynomial (x-r1)(x-r2)(x-r3)(x-r4)=x^4 + a3*x^3 + a2*x^2 + a1*x+ a0
# General formula is called Veita's formula
function vieta_formula(roots)
    a3 = -sum(roots)
    a2 = roots[1]*roots[2] +  roots[1]*roots[3] + roots[1]*roots[4] + roots[2]*roots[3] + roots[2]*roots[4] + roots[3]*roots[4]
    a1 = -(roots[1]*roots[2]*roots[3] + roots[1]*roots[2]*roots[4] + roots[1]*roots[3]*roots[4] + roots[2]*roots[3]*roots[4])
    a0 = prod(roots)
    [a0, a1, a2, a3, 1]
end

@testset "solve general" begin
    r_expected = [1, 2, 3, 4]
    coeffs = complex(vieta_formula(r_expected))

    r_ans = solve_quartic_eq(coeffs)
    r_ans = sort(real.(collect(r_ans)))

    @test r_expected ≈ r_ans
end;

@testset "solve bidratic" begin
    r_quartic = [3.0, 5.3]
    coeffs = complex([prod(r_quartic), 0, -sum(r_quartic), 0, 1])

    r_ans = solve_quartic_eq(coeffs)
    r_ans = sort(real.(collect(r_ans)))

    r_expected = sort([
        -sqrt(r_quartic[1]),
        -sqrt(r_quartic[2]),
        +sqrt(r_quartic[1]),
        +sqrt(r_quartic[2]),
        ])

    @test r_expected ≈ r_ans
end;

@testset "solve zero root" begin
    r_expected = [0.0, 2.1, 3.0, 5.3]
    coeffs = complex(vieta_formula(r_expected))

    r_ans = solve_quartic_eq(coeffs)
    r_ans = sort(real.(collect(r_ans)))

    @test r_expected ≈ r_ans
end;

@testset "solve random" begin
    r_expected = sort(randn(4))
    coeffs = complex(vieta_formula(r_expected))

    r_ans = solve_quartic_eq(coeffs)
    r_ans = sort(real.(collect(r_ans)))

    @test r_expected ≈ r_ans
end;

@testset "solve random complex" begin
    r_expected = randn(4) + im * randn(4)
    sort!(r_expected, by=x->real(x))
    coeffs = vieta_formula(r_expected)

    r_ans = solve_quartic_eq(coeffs)
    r_ans = sort((collect(r_ans)), by=x->real(x))

    @test r_expected ≈ r_ans
end;

