using PolynomialRoots: solve_cubic_eq, solve_quadratic_eq

"""
    solve_quartic_eq([e, d, c, b, a])

Return the roots of `e + dx + cx^2 + bx^3 + ax^4 = 0`.
"""
function solve_quartic_eq(poly::AbstractVector{Complex{T}}) where {T<:AbstractFloat}
    #https://en.wikipedia.org/wiki/Quartic_equation
    @assert length(poly) == 5 "poly must have 5 coefficients"
    if poly[5] == 0
        return solve_cubic_eq(poly[1:4])
    end
    # depressed quartic u^4 + au^2 + bu + c = 0 where x = u - B/4A
    if poly[4] == 0
        c, b, a = poly[[1, 2, 3]] / poly[5]
    else
        E, D, C, B, A = poly
        a = -3B*B/(8A*A) + C/A
        b = B*B*B/(8A*A*A) - B*C/(2A*A) + D/A
        c = -3B*B*B*B/(256*A*A*A*A) + C*B*B/(16*A*A*A) - B*D/(4*A*A) + E/A
    end
    if b == 0 # bidratic
        z1, z2 = solve_quadratic_eq([c, a, 1])
        u = (+sqrt(z1), -sqrt(z1), +sqrt(z2), -sqrt(z2))
    else # general solution
        p = -a*a/12 - c
        q = -a*a*a/108 + a*c/3 - b*b/8
        w = (-q/2 + sqrt(q*q/4 + p*p*p/27))^(1/3)
        y = a/6 + w - p/(3w) 
        z = sqrt(2y - a)
        u = (
            0.5 * (-z + sqrt(-2y - a + 2b/z)),
            0.5 * (-z - sqrt(-2y - a + 2b/z)),
            0.5 * (+z + sqrt(-2y - a - 2b/z)),
            0.5 * (+z - sqrt(-2y - a - 2b/z)),
            )
    end
    x = u .- poly[4]/(4*poly[5])
    x
end

solve_quartic_eq(poly::AbstractVector{<:Number}) = solve_quartic_eq(float.(complex(poly)))

function apply_polynomial(x::Complex{T}, poly::AbstractVector{Complex{T}}) where {T<:AbstractFloat}
    y = 0.0 + 0im
    xp = 1.0 + 0im
    for c in poly
        y += c * xp
        xp *= x
    end
    y
end
