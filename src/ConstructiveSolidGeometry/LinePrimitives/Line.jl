struct Line{T,TP} <: AbstractLinePrimitive{T}
    p1::TP
    p2::TP
    function Line( ::Type{T},
                   p1::Union{CartesianPoint{T}, PlanarPoint{T}},
                   p2::Union{CartesianPoint{T}, PlanarPoint{T}}) where {T}
        new{T,typeof(p1)}(p1, p2)
    end
end

function in(point::Union{PlanarPoint{T}, CartesianPoint{T}}, l::Line{T}) where {T}
    lvec = l.p2 - l.p1
    v = point - l.p1
    isapprox(norm(cross(v, lvec)), 0, atol = geom_atol_zero(T))
end

function distance_to_line(point::Union{PlanarPoint{T}, CartesianPoint{T}},
                          l::Union{Line{T,PlanarPoint{T}}, Line{T,CartesianPoint{T}}}
                          )::T where {T}
    v12 = normalize(l.p2 - l.p1)
    v_point_1 = point - l.p1
    proj_on_v12 = dot(v12,v_point_1)
    return sqrt(abs(dot(v_point_1,v_point_1) - proj_on_v12^2))
end
