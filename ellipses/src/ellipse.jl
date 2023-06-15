include("solve_quartic_eq.jl")

function make_rotation_matrix(θ::T) where {T<:AbstractFloat}
    R = [
        [cos(θ) -sin(θ)]
        [sin(θ) +cos(θ)]
    ]
    R
end

"""
    make_ellipse(a, b, R, [origin]; length=100)
    make_ellipse(a, b, [angle, origin]; length=100)

Return points on an ellipse with semi-axis lengths `a` and `b` rotated by `R` and translated to `origin`.

The points satisfy `(x'/a)^2 + (y'/b)^2 = 1` where `[x'; y'] = R[x + origin[1]; y + origin[2]]`.
"""
function make_ellipse(
    a::T, b::T, R::AbstractMatrix{T}, origin::Vector{T}=zeros(T, 2)
    ; length::Int=1000
    ) where {T<:AbstractFloat}
    t = range(0.0, 2π; length=length)
    x_ellipse = a .* cos.(t)
    y_ellipse = b .* sin.(t)
    X_ellipse = R * transpose([x_ellipse y_ellipse]) .+ origin
    X_ellipse
end

function make_ellipse(
    a::T, b::T, angle::T=0.0, origin::Vector{T}=zeros(T, 2);
    length::Int=1000
    ) where {T<:AbstractFloat}
    R = make_rotation_matrix(angle)
    make_ellipse(a, b, R, origin; length=length)
end

"""
    point_in_ellipse(points, a, b, R, [origin]; length=100)
    point_in_ellipse(points, a, b, [angle, origin]; length=100)

Determine if points lie within or on an ellipse with semi-axis lengths `a` and `b` that is rotated by `R` and translated to `origin`.

True if the point satisfies `(x'/a)^2 + (y'/b)^2 ≤ 1` where `[x'; y'] = inv(R)[x - origin[1]; y - origin[2]]`.
"""
function point_in_ellipse(
    points::AbstractMatrix, a::T, b::T, R::AbstractMatrix{T}, origin::Vector{T}=zeros(T, 2)
    ; atol=1e-6
    ) where {T<:AbstractFloat}
    @assert size(points, 1) == 2 "size(points)=$(size(points)), points must be 2×N"
    points_norm = transpose(R) * (points .- origin)
    inside = (points_norm[1, :] .^ 2 / a^2 + points_norm[2, :] .^ 2 / b^2) .<= (1.0 + atol)
    inside
end

function point_in_ellipse(
    points::AbstractMatrix, a::T, b::T, angle::T=0.0, origin::Vector{T}=zeros(T, 2); options...
    ) where {T<:AbstractFloat}
    R = make_rotation_matrix(angle)
    point_in_ellipse(points, a, b, R, origin; options...)
end

function points_on_ellipse(
    points::AbstractMatrix, a::T, b::T, R::AbstractMatrix{T}, origin::Vector{T}=zeros(T, 2); atol=1e-6
    ) where {T<:AbstractFloat}
    @assert size(points, 1) == 2 "size(points)=$(size(points)), points must be 2×N"
    points_norm = transpose(R) * (points .- origin)
    on_border = isapprox((points_norm[1, :] .^ 2 / a^2 + points_norm[2, :] .^ 2 / b^2), 1; atol=atol)
    on_border
end

function points_on_ellipse(
    points::AbstractMatrix, a::T, b::T, angle::T=0.0, origin::Vector{T}=zeros(T, 2); options...
    ) where {T<:AbstractFloat}
    R = make_rotation_matrix(angle)
    point_in_ellipse(points, a, b, R, origin; options...)
end

"""
    conic_coefficients(a, b, [angle, origin]; length=100)

Return conic section co-efficients for an ellipse.

These satisfy the equation `a0 + a1*x + a2*y + a3*x^2 + a4*xy + a5*y^2 = 0`.
"""
function conic_coefficients(a::Float64, b::Float64, angle::Float64=0.0, origin=[0.0; 0.0])
    A = cos(angle)^2 / a^2 + sin(angle)^2 / b^2
    B = 2 * (1 / a^2 - 1 / b^2) * cos(angle) * sin(angle)
    C = sin(angle)^2 / a^2 + cos(angle)^2 / b^2
    x0 = origin[1]
    y0 = origin[2]

    a0 = A * x0^2 + B * x0 * y0 + C * y0^2 - 1
    a1 = -(2 * A * x0 + B * y0)
    a2 = -(2 * C * y0 + B * x0)
    a3 = A
    a4 = B
    a5 = C
    [a0, a1, a2, a3, a4, a5]
end

"""
    intersect_ellipses(a1, b1, θ1, o1, a2, b2, θ2, o2; atol=1e-6)

The intersection points of 2 overlapping ellipses. This may involve solving a quartic, quadatric or linear equation.
The base simulatenous equations are:
```
    a0 + a1*x + a2*y + a3*x^2 + a4*xy + y^2 = 0
    b0 + b1*x + b2*y + b3*x^2 + b4*xy + y^2 = 0
```

Returns the points and the corresponding number of repeated roots in the polynomials.

If they are the same ellipse then returns nothing.

The absolute tolerance `atol` determines the threshold for significant imaginary components, repeated roots and bounds.
"""
function intersect_ellipses(
    a1::T, b1::T, θ1::T, o1::Vector{T},
    a2::T, b2::T, θ2::T, o2::Vector{T}
    ; atol=1e-6
) where {T<:AbstractFloat}
    coeffs1 = conic_coefficients(a1, b1, θ1, o1)
    coeffs1 /= coeffs1[6]
    coeffs2 = conic_coefficients(a2, b2, θ2, o2)
    coeffs2 /= coeffs2[6]
    d = coeffs1 - coeffs2

    origin_separation = sqrt(sum((o1 - o2).^2))
    # same ellipse or obviously far part
    if all(x -> x == 0, d) || (origin_separation > (max(a1, b1) + max(a2, b2)))
        return Matrix{T}(undef, 2, 0), Int[]
    end

    if d[3] == 0.0 && d[5] == 0.0
        # e.g. both at the origin and concentric
        if d[4] == 0 
            # eg. identical ellipses but horizontally translated
            x_roots = [-d[1] / d[2]]
            n_roots = 1
        else
            x_coeffs = complex(d[[1, 2, 4]])
            x_complex = vcat(solve_quadratic_eq(x_coeffs)...)
            # choose real roots
            almost_real = abs.(imag(x_complex)) .< atol
            x_roots = real.(x_complex[almost_real])
            n_roots = 2
        end
        c0 = coeffs1[1] .+ coeffs1[2] * x_roots + coeffs1[4] * x_roots .* x_roots
        c1 = coeffs1[3] .+ coeffs1[5] * x_roots
        c2 = ones(length(x_roots))
        y_coeffs = complex([c0 c1 c2])
        y1 = solve_quadratic_eq(y_coeffs[1, :])
        if n_roots == 2
            y2 = solve_quadratic_eq(y_coeffs[2, :])
            y_roots = real.(vcat(y1..., y2...))
        else
            y_roots = real.(vcat(y1...))
        end
        x_roots = repeat(x_roots, inner=2)
    else
        f = quartic_coefficients(coeffs1, coeffs2)
        x_complex = vcat(solve_quartic_eq(f)...)
        # choose real roots
        almost_real = abs.(imag(x_complex)) .< atol
        x_roots = real.(x_complex[almost_real])
        y_roots = -(d[1] .+ d[2] * x_roots + d[4] * x_roots .* x_roots) ./ (d[3] .+ d[5] * x_roots)
    end
    points = transpose(hcat(x_roots, y_roots))
    # points in range
    inside = point_in_ellipse(points, a1, b1, θ1, o1; atol=atol)
    points = points[:, inside]
    # remove repeated roots
    unique_idxs, multiplicity = count_repeats(points; atol=atol)
    points[:, unique_idxs], multiplicity
end

"""
    quartic_coefficients(a, b)

Intermediate step for intersecting ellipses. Here `a` and `b` are ellipse co-efficients.

The cofficients of `f0 + f1*x + f2*x^2 + f3*x^3 + f4*x^4 = 0` which solve the simulataneous equations:
```
    a0 + a1*x + a2*y + a3*x^2 + a4*xy + y^2 = 0
    b0 + b1*x + b2*y + b3*x^2 + b4*xy + y^2 = 0
```
Where `(a2 ≠ b2) && (a4 ≠ b4)`.
"""
function quartic_coefficients(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
    @assert length(a) == 6 && length(b) == 6 "vectors must be of length 6"
    a /= a[6]
    b /= b[6]

    d = a - b
    c1 = a[5] * d[3] + a[3] * d[5]
    f0 = a[1] * d[3]^2 - a[3] * d[1] * d[3] + d[1]^2
    f1 = 2 * a[1] * d[3] * d[5] + a[2] * d[3]^2 - a[3] * d[2] * d[3] - c1 * d[1] + 2 * d[1] * d[2]
    f2 = a[1] * d[5]^2 + 2 * a[2] * d[3] * d[5] + a[4] * d[3]^2 - a[3] * d[3] * d[4] - c1 * d[2] - a[5] * d[1] * d[5] + 2 * d[1] * d[4] + d[2]^2
    f3 = a[2] * d[5]^2 + 2 * a[4] * d[3] * d[5] - c1 * d[4] - a[5] * d[2] * d[5] + 2 * d[2] * d[4]
    f4 = a[4] * d[5]^2 - a[5] * d[4] * d[5] + d[4]^2
    [f0, f1, f2, f3, f4]
end

function count_repeats(X::AbstractMatrix{T}; atol=1e-6) where T <: AbstractFloat
    # O(n^2) algorithm. Ill defined → order might affect results
    n = size(X, 2)
    seen = Set{Int}()
    unique_idxs = Int[]
    counts = Int[]
    for i in 1:n
        if i in seen
            continue
        end
        push!(unique_idxs, i)
        m = 1
        for j in (i + 1):n
            if all(abs.(X[:, i] - X[:, j]) .< atol)
                push!(seen, j)
                m += 1
            end
        end
        push!(counts, m)
    end
    unique_idxs, counts
end

"""
    ellipses_intersect_area(a1, b1, θ1, o1, a2, b2, θ2, o2)

The area of intersection of two ellipses. 

There are 5 cases:
1. No intersection if far apart.
2. One ellipse is inside the other.
3. They are identical.
4. Two or three points of intersection, requiring the area of 2 segments.
5. Four points of intersection, requiring 4 segments and 1 quadrilateral.

Reference: "Calculating ellipse overlap areas by Gary B. Hughes and Mohcine Chraibi (2011)" at https://arxiv.org/abs/1106.3787.
"""
function ellipses_intersect_area(
    a1::T, b1::T, θ1::T, o1::Vector{T},
    a2::T, b2::T, θ2::T, o2::Vector{T}
    ; atol=1e-6
) where T <: AbstractFloat

    distance = sqrt(sum((o1 - o2).^2))
    if (distance > (max(a1, b1) + max(a2, b2)))
        return 0.0
    end

    R1 = make_rotation_matrix(θ1)
    R2 = make_rotation_matrix(θ2)
    extrema1 = ellipse_extrema(a1, b1, R1, o1)
    extrema2 = ellipse_extrema(a2, b2, R2, o2)

    p1_in_2 = point_in_ellipse(extrema1, a2, b2, R2, o2)
    p2_in_1 = point_in_ellipse(extrema2, a1, b1, R1, o1)

    if all(p1_in_2)
        return π * a1 * b1
    elseif all(p2_in_1)
        return π * a2 * b2
    end

    points, multiplicity = intersect_ellipses(a1, b1, θ1, o1, a2, b2, θ2, o2; atol=atol)
    num_intersections = length(multiplicity)

    area = 0.0
    if num_intersections == 2 || num_intersections == 3
        if num_intersections == 3 # tangent point is a repeated root
            points = points[:, multiplicity .== 1]
        end

        middle = arc_midpoint(points, a1, b1, R1, o1)
        if point_in_ellipse(middle, a2, b2, R2, o2)[1]
            segment1 = ellipse_segment_area(points, a1, b1, R1, o1)
        else
            segment1 = ellipse_segment_area(reverse(points, dims=2), a1, b1, R1, o1)
        end

        middle = arc_midpoint(points, a2, b2, R2, o2)
        if point_in_ellipse(middle, a1, b1, R1, o1)[1]
            segment2 = ellipse_segment_area(points, a2, b2, R2, o2)
        else
            segment2 = ellipse_segment_area(reverse(points, dims=2), a2, b2, R2, o2)
        end
        area += segment1 + segment2
    elseif num_intersections == 4
        angles = points_to_angles(points, a1, b1, R1, o1)
        idxs = sortperm(angles)
        angles = angles[idxs]
        points = points[:, idxs]
        for i in 1:4
            points_i = points[:, [i, i % 4 + 1]]
            middle = arc_midpoint(points_i, a1, b1, R1, o1)
            if point_in_ellipse(middle, a2, b2, R2, o2)[1]
                area += ellipse_segment_area(points_i, a1, b1, R1, o1)
            else
                area += ellipse_segment_area(points_i, a2, b2, R2, o2)
            end
        end
        points_norm = transpose(R1) * (points .- o1)
        area_quad = 0.5 * abs(
            (points_norm[1, 3] - points_norm[1, 1]) * (points_norm[2, 4] - points_norm[2, 2]) 
            - (points_norm[1, 4] - points_norm[1, 2]) * (points_norm[2, 3] - points_norm[2, 1]) 
        )
        area += area_quad
    end
    area
end

function ellipse_extrema(
    a::T, b::T, R::AbstractMatrix{T}, origin::Vector{T},
) where T <: AbstractFloat
    points = R * [
        -a +a 0.0 0.0;
        0.0 0.0 -b  b
    ] .+ origin
    points
end

function ellipse_extrema(
    a::T, b::T, angle::T, origin::Vector{T},
) where T <: AbstractFloat
    R = make_rotation_matrix(angle)
    ellipse_extrema(a, b, R, origin)
end

"""
    ellipse_segment_area(points, a, b, R, [origin])
    ellipse_segment_area(points, a, b, [θ, origin])

Area of an ellipse segment given by `(t2-t1)ab/2 ± |x1 * y2 - x2 * y1|/2`.
Convention is counter-clockwise so `t2 > t1`.
"""
function ellipse_segment_area(points::Matrix{T}, a::T, b::T, R::AbstractMatrix{T}, origin::Vector{T}=zeros(T, 2)) where T <: AbstractFloat
    # algorithm 2.4 from Calculating ellipse overlap areas by Gary B. Hughes and Mohcine Chraibi
    @assert size(points) == (2, 2) "size(points)=$(size(points)), points must be 2×2"

    t = points_to_angles(points, a, b, R, origin)
    t[1] = t[1] <= t[2] ? t[1] : t[1] - 2π # require t1 < t2 for integral

    area_sector = 0.5 * (t[2] - t[1]) * a * b
    points_norm = transpose(R) * (points .- origin)
    x1 = points_norm[1, 1]
    y1 = points_norm[2, 1]
    x2 = points_norm[1, 2]
    y2 = points_norm[2, 2]
    area_triangle = 0.5 * abs(x1 * y2 - x2 * y1)

    area_sector + sign(t[2] - t[1] - π) * area_triangle
end

function ellipse_segment_area(points::Matrix{T}, a::T, b::T, angle::T=0.0, origin::Vector{T}=zeros(T, 2)) where T <: AbstractFloat
    R = make_rotation_matrix(angle)
    ellipse_segment_area(points, a, b, R, origin)
end

"""
    points_to_angles(points, a, b, R, origin)
    points_to_angles(points, a, b, [θ, origin])

Parametric angles of an ellipse for the equation:
```
[x; y] = R * [a*cos(t); b*sin(t)] + origin
```
"""
function points_to_angles(points::Matrix{T}, a::T, b::T, R::AbstractMatrix{T}, origin::Vector{T}=zeros(T, 2)) where T
    points_norm = transpose(R) * (points .- origin)
    x = clamp.(points_norm[1, :] / a, -1, 1)
    t = acos.(x) 
    t = ifelse.(points_norm[2, :] .>= 0, t, 2π .- t)
    t
end

function points_to_angles(points::Matrix{T}, a::T, b::T, angle::T=0.0, origin::Vector{T}=zeros(T, 2)) where T
    R = make_rotation_matrix(angle)
    points_to_angles(points, a, b, R, origin)
end

"""
    arc_midpoint(points, a, b, R, origin)
    arc_midpoint(points, a, b, [θ, origin])

The midpoint of the arc of the ellipse between two points going in a counter-clockwise direction.
"""
function arc_midpoint(points::Matrix{T}, a::T, b::T, R::AbstractMatrix{T}, origin::Vector{T}=zeros(T, 2)) where T
    @assert size(points) == (2, 2) "size(points)=$(size(points)), points must be 2×2"
    t = points_to_angles(points, a, b, R, origin)
    t[1] = t[1] <= t[2] ? t[1] : t[1] - 2π # require θ1 < θ2 for integral
    t_mid = (t[1] + t[2]) / 2
    p_mid = R * reshape([a * cos(t_mid); b * sin(t_mid)], :, 1) + origin
    p_mid
end

function arc_midpoint(points::Matrix{T}, a::T, b::T, angle::T=0.0, origin::Vector{T}=zeros(T, 2)) where T
    R = make_rotation_matrix(angle)
    arc_midpoint(points, a, b, R, origin)
end
