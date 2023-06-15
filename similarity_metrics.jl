# Polygons
using PolygonAlgorithms: matrix_to_points
using PolygonAlgorithms: convex_hull, intersect_convex, area_polygon, GrahamScanAlg

# Ellipse
include("ellipses/src/ellipse.jl")
include("ellipses/src/eigenvalues.jl")

# Hungarian
using Hungarian

function similarities_polygons(points1::AbstractVector, points2::AbstractVector)   
    hull_idxs1 = convex_hull(points1);
    hull1 = points1[hull_idxs1]

    hull_idxs2 = convex_hull(points2);
    hull2 = points2[hull_idxs2]
        
    intersection_points = intersect_convex(hull1, hull2);
    area1 = area_polygon(hull1)
    area_inter = area_polygon(intersection_points)
    area2 = area_polygon(hull2)
    IoU = area_inter / (area1 + area2 - area_inter)
    IoU
end

function similarities_polygons(points1::AbstractMatrix, points2::AbstractMatrix)
    points1_vec = matrix_to_points(points1)
    points2_vec = matrix_to_points(points2)
    similarities_polygons(points1_vec, points2_vec)
end

function similarities_ellipses(points1::AbstractMatrix, points2::AbstractMatrix; scale=2.45)
    C1 = cov(points1, dims=2)
    C2 = cov(points2, dims=2)
    
    λ1s, W1 = calc_eigvals_2x2(C1)
    λ2s, W2 = calc_eigvals_2x2(C2)
    
    a1 = scale * sqrt(λ1s[1])
    b1 = scale * sqrt(λ1s[2])
    angle1 = atan(W1[2, 1], W1[1, 1])
    origin1 = vec(mean(points1, dims=2))

    a2 = scale * sqrt(λ2s[1])
    b2 = scale * sqrt(λ2s[2])
    angle2 = atan(W2[2, 1], W2[1, 1])
    origin2 = vec(mean(points2, dims=2))
    
    area_inter = ellipses_intersect_area(a1, b1, angle1, origin1, a2, b2, angle2, origin2)
    area1 = π * a1 * b1
    area2 = π * a2 * b2
    IoU = area_inter / (area1 + area2 - area_inter)
    IoU
end

function make_distance_matrix_vectorised(points1::AbstractMatrix, points2::AbstractMatrix)
    square1 = transpose(sum(points1 .* points1, dims=1))
    square2 = sum(points2 .* points2, dims=1)
    distance_matrix = square1 .- 2 * transpose(points1) * points2 .+ square2
    distance_matrix
end

function similarities_hungarian(points1, points2)
    weights = make_distance_matrix_vectorised(points1, points2)
    assignments, cost = hungarian(weights)
    cost / min(size(weights)...)
end