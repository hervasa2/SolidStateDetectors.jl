PlanarSurface{T} = Union{CylindricalAnnulus{T}, ToroidalAnnulus{T}, ConalPlane{T}, AbstractConstructiveSurface{T}}

struct PlanarSurfaceDifference{T, A, B} <: AbstractConstructiveSurface{T}
    a::A
    b::B
    lines::Array{AbstractPrimitive}
    function PlanarSurfaceDifference( ::Type{T}, a::PlanarSurface{T}, b::PlanarSurface{T}) where {T}
        empty = new{T,typeof(a),typeof(b)}(a, b, AbstractPrimitive[])
        #lines = get_decomposed_surfaces(empty)
        lines = AbstractPrimitive[]
        new{T,typeof(a),typeof(b)}(a, b, lines)
    end
end

(-)(a::A, b::B) where {T, A <: PlanarSurface{T}, B <: PlanarSurface{T}} = PlanarSurfaceDifference{T,A,B}(a, b)

in(p::AbstractCoordinatePoint, csg::PlanarSurfaceDifference) = in(p, csg.a) && !in(p, csg.b)

#not for unions
sample(ps::AbstractConstructiveSurface{T}, step::Union{Real, NTuple{3,Int}}) where {T} = filter!(x -> x in ps, sample(ps.a, step))

Plane(ps::AbstractConstructiveSurface) = Plane(ps.a)

function distance_to_surface(point::AbstractCoordinatePoint{T}, ps::AbstractConstructiveSurface{T})::T where {T}
    pct = CartesianPoint(point)
    plane = Plane(ps)
    projection, distance_to_plane = project_on_plane(pct, plane)
    if projection in ps
        return distance_to_plane
    else
        d = T(Inf)
        for line in ps.lines
            dist = distance_to_line(point,line)
            dist < d ? d = dist : nothing
        end
        return NaN#hypot(d, distance_to_plane)
    end
end
