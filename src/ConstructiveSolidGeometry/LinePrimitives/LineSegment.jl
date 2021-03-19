struct LineSegment{T,TP} <: AbstractLinePrimitive{T}
    p1::TP
    p2::TP
    function LineSegment( ::Type{T},
                   p1::Union{CartesianPoint{T}, PlanarPoint{T}},
                   p2::Union{CartesianPoint{T}, PlanarPoint{T}}) where {T}
        new{T,typeof(p1)}(p1, p2)
    end
end

Line(l::LineSegment{T}) where {T} = Line(T,l.p1,l.p2)
LineSegment(l::Line{T}) where {T} = LineSegment(T,l.p1,l.p2)

translate(l::LineSegment{T,PlanarPoint{T}}, t::PlanarVector{T}) where {T} = LineSegment(T, l.p1 + t, l.p2 + t)

get_line_vector(l::LineSegment{T,PlanarPoint{T}}) where {T} = normalize(PlanarVector{T}(l.p2 - l.p1))
get_line_vector(l::LineSegment{T,CartesianPoint{T}}) where {T} = normalize(CartesianPointVector{T}(l.p2 - l.p1))

function get_surface_vector(l::LineSegment{T,PlanarPoint{T}}) where {T}
    v = get_line_vector(l)
    return PlanarVector{T}(-v.v, v.u)
end

get_midpoint(l::LineSegment{T}) where {T} = (l.p2 + l.p1)/2

function in(point::Union{PlanarPoint{T}, CartesianPoint{T}}, l::LineSegment{T}) where {T}
    if point in (l.p2, l.p1)
        return true
    else
        lv = l.p2 - l.p1
        length = norm(lv)
        lvec = (lv)/length
        v = point - l.p1
        tol = geom_atol_zero(T)
        return isapprox(norm(cross(v, lvec)), T(0), atol = tol) && -tol ≤ dot(v, lvec) ≤ length + tol
    end
end

function distance_to_line(point::Union{PlanarPoint{T}, CartesianPoint{T}},
                          l::Union{LineSegment{T,PlanarPoint{T}}, LineSegment{T,CartesianPoint{T}}}
                          )::T where {T}
    v12 = get_line_vector(l)
    v_point_1 = point - l.p1
    proj_on_v12 = dot(v12,v_point_1)
    if proj_on_v12 ≤ T(0)
        return norm(l.p1 - point)
    else
        v_point_2 = point - l.p2
        if dot(v12,v_point_2) ≥ T(0)
            return norm(l.p2 - point)
        else
            return sqrt(abs(dot(v_point_1,v_point_1) - proj_on_v12^2))
        end
    end
end

function get_Ax_By_C_line_representation(l::Union{Line{T,PlanarPoint{T}}, LineSegment{T,PlanarPoint{T}}}) where {T}
    # Ax + By = C where C = 0 or C = 1, points given are the same will assume one is (0,0)
    M = transpose(hcat(l.p1,l.p2))
    if det(M) == 0 #non trivial solution to homogeneus eq Mx = 0
        if l.p1.u == 0
            return l.p2.u == 0 ? (T(1), T(0), T(0)) : (-l.p2.v/l.p2.u, T(1), T(0))
        else
            return -l.p1.v/l.p1.u, T(1), T(0)
        end
    else
        if l.p1.u == l.p2.u && l.p1.v ≠ l.p2.v
            return T(1), T(0), l.p1.u
        else
            b = PlanarPoint{T}(1,1) #solving Mx = b
            x = inv(M)*b
            return x[1], x[2], T(1)
        end
    end
end

function cut(l::LineSegment{T}, point::Union{CartesianPoint{T}, PlanarPoint{T}}, ::Val{:point}) where {T}
    if point == l.p1 || point == l.p2
        return LineSegment{T}[l]
    else
        if point in l
            return LineSegment{T}[
                            LineSegment(T, l.p1, point),
                            LineSegment(T, point, l.p2)
                         ]
        else
            return LineSegment{T}[l]
        end
    end
end
