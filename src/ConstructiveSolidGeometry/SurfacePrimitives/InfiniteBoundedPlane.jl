#only meant to be used as second field of PlanarSurfaceDifference and Intersection
struct InfiniteBoundedPlane{T} <: AbstractPlanarSurfacePrimitive{T}
    plane::Plane{T} #note that plane should have orthogonal spanning vectors or given line will be in different coordinate system
    line::Line{T,PlanarPoint{T}} #boundary
    dir::PlanarVector{T} #in which direction does the infinite plane span to
end

Plane(ip::InfiniteBoundedPlane{T}) where {T} = ip.plane

in(point::PlanarPoint, ip::InfiniteBoundedPlane) = dot(point - ip.line.p1, ip.dir) ≥ 0

function in(point::AbstractCoordinatePoint, ip::InfiniteBoundedPlane)
    pct = CartesianPoint(point)
    dot(get_planar_point(pct, ip.plane) - ip.line.p1, ip.dir) ≥ 0 && on_infinite_plane(pct, ip.plane)
end

get_decomposed_lines(ip::InfiniteBoundedPlane) = AbstractLinePrimitive[ip.line]
