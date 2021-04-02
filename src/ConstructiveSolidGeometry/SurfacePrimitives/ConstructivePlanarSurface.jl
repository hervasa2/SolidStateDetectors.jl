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


function attempt_reduce_to_primitive(ps::Union{PlanarSurfaceDifference{T, <:Any, <:Any, Val{:φ}}, PlanarSurfaceIntersection{T, <:Any, <:Any, Val{:φ}}}) where {T}
    φ = get_plane_φ(ps)
    Arcs = filter(l -> typeof(l) <: Arc{T}, ps.lines)
    Segs = filter(l -> typeof(l) <: Line{T, <:Any, Val{:seg}}, ps.lines)
    nSegs = length(Segs)
    nArcs = length(Arcs)
    if nArcs == 0
        vertices = unique(append!([geom_round(l.p1) for l in ps.lines], [geom_round(l.p2) for l in ps.lines]))
        sort!(vertices, by = p -> p.v)
        if length(vertices) == 3
            zMin, zMax = vertices[1].v, vertices[3].v
            if vertices[1].v == vertices[2].v
                rtopMin, rtopMax = vertices[3].u, vertices[3].u
                rbotMin, rbotMax = minmax(vertices[1].u, vertices[2].u)
                return ConalPlane(rbotMin, rbotMax, rtopMin, rtopMax, φ, zMin, zMax)
            elseif vertices[2].v == vertices[3].v
                rtopMin, rtopMax = minmax(vertices[2].u, vertices[3].u)
                rbotMin, rbotMax = vertices[1].u, vertices[1].u
                return ConalPlane(rbotMin, rbotMax, rtopMin, rtopMax, φ, zMin, zMax)
            else
                return ps
            end
        elseif length(vertices) == 4
            zMin, zMax = vertices[1].v, vertices[4].v
            if vertices[1].v == vertices[2].v && vertices[3].v == vertices[4].v
                rtopMin, rtopMax = minmax(vertices[3].u, vertices[4].u)
                rbotMin, rbotMax = minmax(vertices[1].u, vertices[2].u)
                return ConalPlane(rbotMin, rbotMax, rtopMin, rtopMax, φ, zMin, zMax)
            else
                return ps
            end
        else
            return ps
        end
    elseif nArcs == 1 && nSegs == 0
        return ps
    elseif nArcs == 1 && nSegs == 2
        return ps
    elseif nArcs == 2 && nSegs == 2
        return ps
    else
        return ps
    end
end

#not for unions
sample(ps::AbstractConstructivePlanarSurface{T}, step::Union{Real, NTuple{3,Int}}) where {T} = filter!(x -> x in ps, sample(ps.a, step))

function sample_border(ps::AbstractConstructivePlanarSurface{T}, sampling) where {T}
    samples = [ point
    for line in ps.lines
    for point in sample(line, sampling)  ]
end

perimeter(ps::AbstractConstructivePlanarSurface{T}) where {T} = sum([ perimeter(line) for line in ps.lines ])

function get_nodes(ps::AbstractConstructivePlanarSurface{T}, n_arc::Real) where {T}
    Arcs = filter(l -> typeof(l) <: Arc{T}, ps.lines)
    Segs = filter(l -> typeof(l) <: Line{T, <:Any, Val{:seg}}, ps.lines)
    nSegs = length(Segs)
    nArcs = length(Arcs)
    n_seg = nArcs == 0 ? n_arc : max(2, 1 + Int(floor(n_arc*nArcs/nSegs)))
    arc_nodes = unique([ geom_round(point)
                        for arc in Arcs
                        for point in get_nodes(arc, n_arc) ]
                      )
    seg_nodes = unique([ geom_round(point)
                        for seg in Segs
                        for point in get_nodes(seg, n_seg) ]
                      )
    unique(append!(seg_nodes, arc_nodes))
end

function get_midpoint(ps::Union{PlanarSurfaceDifference{T, <:Any, <:Any, Val{:φ}}, PlanarSurfaceIntersection{T, <:Any, <:Any, Val{:φ}}}) where {T}
    if length(ps.lines) ≥ 1
        φ = get_plane_φ(ps)
        line = ps.lines[1]
        tol = 10*geom_atol_zero(T)
        pt = get_midpoint(line)
        n = get_surface_vector(line, pt)
        pt_pos = pt + tol*n
        pt_neg = pt - tol*n
        return pt_pos in ps ? CartesianPoint(CylindricalPoint{T}(pt_pos.u, φ, pt_pos.v)) :
            CartesianPoint(CylindricalPoint{T}(pt_neg.u, φ, pt_neg.v))
    else
        @error "Constructive surface has no lines!"
    end
end

Plane(ps::AbstractConstructivePlanarSurface) = Plane(ps.a)

get_surface_vector(ps::AbstractConstructivePlanarSurface, point::AbstractCoordinatePoint) = get_surface_vector(ps.a, point)

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

function merge(ps1::AbstractConstructivePlanarSurface, ps2::AbstractConstructivePlanarSurface)
    if ps1 == ps2
        return ps1, true
    elseif get_plane_φ(ps1) == get_plane_φ(ps2) && length(ps1.lines) == length(ps2.lines) == sum([l in ps2.lines for l in ps1.lines])
        return ps1, true
    else
        return ps1, false
    end
end
