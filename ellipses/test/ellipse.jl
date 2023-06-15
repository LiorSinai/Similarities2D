using Test

function coeff_error_test(a::Vector{T}, x::T, y::T) where T 
    errors = a[1] + a[2] * x + a[3] * y + a[4] * x * x + a[5] * x * y + a[6] * ye * ye
    maximum(abs.(errors)) < 1e-8
end

function coeff_error_test(a::Vector{T}, x::Vector{T}, y::Vector{T}) where T 
    errors = a[1] .+ a[2] .* x + a[3] .* y + a[4] * x .* x + a[5] * x .* y + a[6] * y .* y
    maximum(abs.(errors)) < 1e-8
end

function validate_points(
    points::AbstractMatrix,
    a1::T, b1::T, θ1::T, o1::Vector{T},
    a2::T, b2::T, θ2::T, o2::Vector{T} 
    ) where {T<:AbstractFloat}
    points_on_1 = points_on_ellipse(points, a1, b1, θ1, o1);
    points_on_2 = points_on_ellipse(points, a2, b2, θ2, o2);
    all(points_on_1) && all(points_on_2)
end

@testset "Coefficients" verbose=true begin
    @testset "basic" begin
        r1 = 2.0
        r2 = 3.0
        X = make_ellipse(r1, r2, length=1000)
        x = X[1, :]
        y = X[2, :]

        coeffs = conic_coefficients(r1, r2)
        @test coeff_error_test(coeffs, x, y)
    end;

    @testset "circle" begin
        r1 = 2.0
        r2 = 2.0
        X = make_ellipse(r1, r2, length=1000)
        x = X[1, :]
        y = X[2, :]

        coeffs = conic_coefficients(r1, r2)
        @test coeff_error_test(coeffs, x, y)
    end;

    @testset "rotated" begin
        r1 = 2.0
        r2 = 3.0
        θ = 30 * π / 180
        X = make_ellipse(r1, r2, θ, length=1000)
        x = X[1, :]
        y = X[2, :]

        coeffs = conic_coefficients(r1, r2, θ)
        @test coeff_error_test(coeffs, x, y)
    end;

    @testset "translated" begin
        r1 = 2.0
        r2 = 3.0
        θ = 0.0
        t = [-1.5, 3.2]
        X = make_ellipse(r1, r2, θ, t, length=1000)
        x = X[1, :]
        y = X[2, :]

        coeffs = conic_coefficients(r1, r2, θ, t)
        @test coeff_error_test(coeffs, x, y)
    end;

    @testset "general" begin
        r1 = 2.0
        r2 = 3.0
        θ = 30  * π / 180
        t = [-1.5, 3.2]
        X = make_ellipse(r1, r2, θ, t, length=1000)
        x = X[1, :]
        y = X[2, :]

        coeffs = conic_coefficients(r1, r2, θ, t)
        @test coeff_error_test(coeffs, x, y)
    end;
end

@testset "intersection points" verbose=true begin
    params1 = (1.0, 1.5, 30 * π / 180, [0.0; 0.0])
    params2 = (2.0, 0.5, 30 * π / 180)

    @testset "same" begin 
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2, origin2 = params1
        intersection_points, multiplicity = intersect_ellipses(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        @test isempty(intersection_points)
    end

    @testset "zero" begin
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2 = params2
        origin2 = [3.0; 0.0]
        intersection_points, multiplicity = intersect_ellipses(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        
        @test isempty(intersection_points)
    end

    @testset "one" begin
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2 = params2
        origin2 = 3 * [cos(θ2); sin(θ2)] 
        intersection_points, multiplicity = intersect_ellipses(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        
        expected_points = [0.8660254037844386; 0.5 ;;]
        expected_multiplicity = [2]
        @test all(intersection_points .≈ expected_points)
        @test all(multiplicity .== expected_multiplicity)
        @test validate_points(intersection_points, a1, b1, θ1, origin1, a2, b2, θ2, origin2)
    end

    @testset "two" begin
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2 = params2
        origin2 = [1.0; 1.0]
        intersection_points, multiplicity = intersect_ellipses(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        
        expected_points = [
            0.9256785429537308  0.29089511609154106;  
            0.38968029397483067 1.1464339992372108
        ]
        expected_multiplicity = [1, 1]
        @test all(intersection_points .≈ expected_points)
        @test all(multiplicity .== expected_multiplicity)
        @test validate_points(intersection_points, a1, b1, θ1, origin1, a2, b2, θ2, origin2)
    end

    @testset "three" begin
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2 = params2
        origin2 = [cos(θ2); sin(θ2)] 
        intersection_points, multiplicity = intersect_ellipses(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        
        expected_points = [
            -0.8660254037844387 1.066436176204412 0.5666402995033866
            -0.5000000000000001  0.03859264549876838 0.9042644973583747
        ]
        expected_multiplicity = [2, 1, 1]
        @test all(intersection_points .≈ expected_points)
        @test all(multiplicity .== expected_multiplicity)
        @test validate_points(intersection_points, a1, b1, θ1, origin1, a2, b2, θ2, origin2)
    end

    @testset "four" begin
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2 = params2
        origin2 = [0.0; 0.0] 
        intersection_points, multiplicity = intersect_ellipses(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        
        expected_points = [
            1.047656187624245 -1.047656187624245 0.6085011547974053 -0.6085011547974053
            0.09777202910592502 -0.09777202910592502 0.8584108583615897 -0.8584108583615897
        ]
        expected_multiplicity = [1, 1, 1, 1]
        @test all(intersection_points .≈ expected_points)
        @test all(multiplicity .== expected_multiplicity)
        @test validate_points(intersection_points, a1, b1, θ1, origin1, a2, b2, θ2, origin2)
    end

    @testset "concentric" begin
        a1, b1, θ1, origin1 = 3.0, 2.0, 0.0, [0.0, 0.0]
        a2, b2, θ2, origin2 = 2.0, 4.0, 0.0, [0.0, 0.0] 

        intersection_points, multiplicity = intersect_ellipses(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        
        x = 1.8371173070873836
        y = 1.5811388300841898
        expected_points = [
            +x +x -x -x
            -y +y -y +y
        ]
        expected_multiplicity = [1, 1, 1, 1]
        @test all(intersection_points .≈ expected_points)
        @test all(multiplicity .== expected_multiplicity)
        @test validate_points(intersection_points, a1, b1, θ1, origin1, a2, b2, θ2, origin2)
    end

    @testset "horizontally translated" begin
        a1, b1, θ1, origin1 = 3.0, 2.0, 0.0, [0.0, 0.0]
        a2, b2, θ2, origin2 = a1, b1, θ1, origin1 + [2; 0]

        intersection_points, multiplicity = intersect_ellipses(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        
        x = 1
        y = 1.8856180831641267
        expected_points = [
            +x +x
            -y +y
        ]
        expected_multiplicity = [1, 1]
        @test all(intersection_points .≈ expected_points)
        @test all(multiplicity .== expected_multiplicity)
        @test validate_points(intersection_points, a1, b1, θ1, origin1, a2, b2, θ2, origin2)
    end

    @testset "roots outside bounds" begin
        a1, b1, θ1, origin1 = 2.5, 1.5, 0.0, [0.0, 0.0]
        a2, b2, θ2, origin2 = 2.0, 1.0, 0.0, [0.0, 0.0] 

        intersection_points, multiplicity = intersect_ellipses(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        @test isempty(intersection_points)
        @test validate_points(intersection_points, a1, b1, θ1, origin1, a2, b2, θ2, origin2)
    end
end

@testset "intersection area" verbose=true begin
    params1 = (1.0, 1.5, 30 * π / 180, [0.0; 0.0])
    params2 = (2.0, 0.5, 30 * π / 180)

    @testset "same" begin 
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2, origin2 = params1
        area = ellipses_intersect_area(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        expected = π * a1 * b1
        @test area ≈ expected
    end

    @testset "inside" begin 
        a1, b1, θ1, origin1 = 2.5, 1.5, 0.0, [0.0, 0.0]
        a2, b2, θ2, origin2 = 2.0, 1.0, 0.0, [0.0, 0.0] 
        area = ellipses_intersect_area(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        expected = π * a2 * b2
        @test area ≈ expected
    end

    @testset "zero" begin
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2 = params2
        origin2 = [3.0; 0.0]
        area = ellipses_intersect_area(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        @test area ≈ 0.0
    end

    @testset "one" begin
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2 = params2
        origin2 = 3 * [cos(θ2); sin(θ2)] 
        area = ellipses_intersect_area(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        @test area ≈ 0.0
    end

    @testset "two" begin
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2 = params2
        origin2 = [1.0; 1.0]
        area = ellipses_intersect_area(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        
        expected = 1.1582030253873428
        @test area ≈ expected
    end

    @testset "three" begin
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2 = params2
        origin2 = [cos(θ2); sin(θ2)] 
        area = ellipses_intersect_area(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        
        expected = 1.5519639970755992
        @test area ≈ expected
    end

    @testset "four" begin
        a1, b1, θ1, origin1 = params1
        a2, b2, θ2 = params2
        origin2 = [0.0; 0.0] 
        area = ellipses_intersect_area(a1, b1, θ1, origin1, a2, b2, θ2, origin2)
        
        expected = 1.8883284233161728
        @test area ≈ expected
    end
end