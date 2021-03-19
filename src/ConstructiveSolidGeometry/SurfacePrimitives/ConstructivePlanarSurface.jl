struct PlanarSurfaceDifference{T, A, B, TP} <: AbstractConstructivePlanarSurface{T}
    a::A
    b::B
    planetype::TP
    lines::Array{AbstractLinePrimitive}
    function PlanarSurfaceDifference( ::Type{T},
                                     a::Union{AbstractPlanarSurfacePrimitive{T}, AbstractConstructivePlanarSurface{T}},
                                     b::Union{AbstractPlanarSurfacePrimitive{T}, AbstractConstructivePlanarSurface{T}},
                                     planetype::Union{Val{:φ}, Val{:z}, Val{:arbitrary}}) where {T}
        empty = new{T,typeof(a),typeof(b),typeof(planetype)}(a, b, planetype, AbstractLinePrimitive[])
        lines = get_decomposed_lines(empty)
        new{T,typeof(a),typeof(b),typeof(planetype)}(a, b, planetype, lines)
    end
end

(-)(a::A, b::B) where {T,
                            A <: Union{AbstractPlanarSurfacePrimitive{T}, AbstractConstructivePlanarSurface{T}},
                            B <: Union{AbstractPlanarSurfacePrimitive{T}, AbstractConstructivePlanarSurface{T}}
                       } = PlanarSurfaceDifference{T,A,B}(a, b)

in(p::Union{AbstractCoordinatePoint, AbstractPlanarPoint}, csg::PlanarSurfaceDifference) = in(p, csg.a) && !in(p, csg.b)

struct PlanarSurfaceIntersection{T, A, B, TP} <: AbstractConstructivePlanarSurface{T}
    a::A
    b::B
    planetype::TP
    lines::Array{AbstractLinePrimitive}
    function PlanarSurfaceIntersection( ::Type{T},
                                     a::Union{AbstractPlanarSurfacePrimitive{T}, AbstractConstructivePlanarSurface{T}},
                                     b::Union{AbstractPlanarSurfacePrimitive{T}, AbstractConstructivePlanarSurface{T}},
                                     planetype::Union{Val{:φ}, Val{:z}, Val{:arbitrary}}) where {T}
        empty = new{T,typeof(a),typeof(b),typeof(planetype)}(a, b, planetype, AbstractLinePrimitive[])
        lines = get_decomposed_lines(empty)
        new{T,typeof(a),typeof(b),typeof(planetype)}(a, b, planetype, lines)
    end
end

(&)(a::A, b::B) where {T,
                            A <: Union{AbstractPlanarSurfacePrimitive{T}, AbstractConstructivePlanarSurface{T}},
                            B <: Union{AbstractPlanarSurfacePrimitive{T}, AbstractConstructivePlanarSurface{T}}
                       } = PlanarSurfaceIntersection{T,A,B}(a, b)

in(p::Union{AbstractCoordinatePoint, AbstractPlanarPoint}, csg::PlanarSurfaceIntersection) = in(p, csg.a) && in(p, csg.b)

get_plane_φ(ps::Union{PlanarSurfaceDifference{<:Any, <:Any, <:Any, Val{:φ}}, PlanarSurfaceIntersection{<:Any, <:Any, <:Any, Val{:φ}}}) = get_plane_φ(ps.a)

#not for unions
sample(ps::AbstractConstructivePlanarSurface{T}, step::Union{Real, NTuple{3,Int}}) where {T} = filter!(x -> x in ps, sample(ps.a, step))

Plane(ps::AbstractConstructivePlanarSurface) = Plane(ps.a)

get_surface_vector(ps::AbstractConstructivePlanarSurface) = get_surface_vector(ps.a)

function distance_to_surface(point::PlanarPoint{T}, ps::AbstractConstructivePlanarSurface{T})::T where {T}
    if point in ps
        return T(0)
    else
        d = T(Inf)
        for line in ps.lines
            dist = distance_to_line(point,line)
            dist < d ? d = dist : nothing
        end
        return d
    end
end

function distance_to_surface(point::AbstractCoordinatePoint{T}, ps::AbstractConstructivePlanarSurface{T})::T where {T}
    pct = CartesianPoint(point)
    plane = Plane(ps)
    projection, distance_to_plane = project_on_plane(pct, plane)
    planarpoint = get_planar_point(projection + plane.p1, plane)
    if planarpoint in ps
        return distance_to_plane
    else
        d = T(Inf)
        for line in ps.lines
            dist = distance_to_line(planarpoint,line)
            dist < d ? d = dist : nothing
        end
        return hypot(d, distance_to_plane)
    end
end

PlanarSurfaceAtφ{T} = Union{
                                ToroidalAnnulus{T},
                                ConalPlane{T},
                                PlanarSurfaceDifference{T, <:Any, <:Any, Val{:z}},
                                PlanarSurfaceIntersection{T, <:Any, <:Any, Val{:z}}
                            }

PlanarSurfaceAtz{T} = Union{
                                CylindricalAnnulus{T},
                                PlanarSurfaceDifference{T, <:Any, <:Any, Val{:z}},
                                PlanarSurfaceIntersection{T, <:Any, <:Any, Val{:z}}
                            }
