struct Line{T,TP,TL} <: AbstractLinePrimitive{T}
    p1::TP
    p2::TP
    linetype::TL
    function Line( ::Type{T},
                   p1::Union{CartesianPoint{T}, PlanarPoint{T}},
                   p2::Union{CartesianPoint{T}, PlanarPoint{T}},
                   linetype::Union{Val{:inf}, Val{:ray}, Val{:seg}}) where {T}
        new{T,typeof(p1),typeof(linetype)}(p1, p2)
    end
end

function Line( ::Type{T},
               p1::Union{CartesianPoint{T}, PlanarPoint{T}},
               p2::Union{CartesianPoint{T}, PlanarPoint{T}}) where {T} #if (inf)Line there are no endpoints
    Line(T, p1, p2, Val(:inf))
end

function Ray( ::Type{T},
               p1::Union{CartesianPoint{T}, PlanarPoint{T}},
               p2::Union{CartesianPoint{T}, PlanarPoint{T}}) where {T} #if Ray endpoint is p1
    Line(T, p1, p2, Val(:ray))
end

function LineSegment( ::Type{T},
               p1::Union{CartesianPoint{T}, PlanarPoint{T}},
               p2::Union{CartesianPoint{T}, PlanarPoint{T}}) where {T} #if LineSegment endpoints are p1, and p2
    Line(T, p1, p2, Val(:seg))
end

Line(l::Line{T}) where {T} = Line(T, l.p1, l.p2)

get_line_vector(l::Line{T,PlanarPoint{T}}) where {T} = normalize(PlanarVector{T}(l.p2 - l.p1))
get_line_vector(l::Line{T,CartesianPoint{T}}) where {T} = normalize(CartesianVector{T}(l.p2 - l.p1))

function in(point::Union{PlanarPoint, CartesianPoint}, l::Line{T, <:Any, Val{:inf}}) where {T}
    lvec = l.p2 - l.p1
    v = point - l.p1
    isapprox(norm(cross(v, lvec)), 0, atol = geom_atol_zero(T))
end

function in(point::Union{PlanarPoint, CartesianPoint}, l::Line{T, <:Any, Val{:ray}}) where {T}
    if point == l.p1
        return true
    else
        lvec = l.p2 - l.p1
        v = point - l.p1
        tol = geom_atol_zero(T)
        return isapprox(norm(cross(v, lvec)), T(0), atol = tol) && -tol ≤ dot(v, lvec)
    end
end

function in(point::Union{PlanarPoint, CartesianPoint}, l::Line{T, <:Any, Val{:seg}}) where {T}
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

get_midpoint(l::Line{T}) where {T} = (l.p2 + l.p1)/2

translate(l::Line{T}, t::Union{CartesianVector, PlanarVector}) where {T} = Line(T, l.p1 + t, l.p2 + t, l.linetype)

function get_surface_vector(l::Line{T,PlanarPoint{T}}) where {T}
    v = get_line_vector(l)
    return PlanarVector{T}(-v.v, v.u)
end

function distance_to_line(point::Union{PlanarPoint{T}, CartesianPoint{T}}, l::Line{T, <:Any, Val{:inf}})::T where {T}
    v12 = normalize(l.p2 - l.p1)
    v_point_1 = point - l.p1
    proj_on_v12 = dot(v12,v_point_1)
    return sqrt(abs(dot(v_point_1,v_point_1) - proj_on_v12^2))
end

function distance_to_line(point::Union{PlanarPoint{T}, CartesianPoint{T}}, l::Line{T, <:Any, Val{:ray}})::T where {T}
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

function distance_to_line(point::Union{PlanarPoint{T}, CartesianPoint{T}}, l::Line{T, <:Any, Val{:seg}})::T where {T}
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

function get_Ax_By_C_line_representation(l::Line{T,PlanarPoint{T}}) where {T}
    # Ax + By = C where C = 0 or C = 1, points given are the same will assume one is (0,0)
    M = transpose(hcat(l.p1,l.p2))
    if det(M) == 0 #non trivial solution to homogeneus eq Mx = 0
        if l.p1.u == 0
            return l.p2.u == 0 ? (T(1), T(0), T(0)) : (-l.p2.v/l.p2.u, T(1), T(0))
        else
            return -l.p1.v/l.p1.u, T(1), T(0)
        end
    else
        if isapprox(l.p1.u, l.p2.u, atol = geom_atol_zero(T)) && l.p1.v ≠ l.p2.v
            return T(1), T(0), l.p1.u
        elseif isapprox(l.p1.v, l.p2.v, atol = geom_atol_zero(T)) && l.p1.u ≠ l.p2.u
            return T(0), T(1), l.p1.v
        else
            b = PlanarPoint{T}(1,1) #solving Mx = b
            x = inv(M)*b
            return x[1], x[2], T(1)
        end
    end
end

function cut(l::Line{T,<:Any,Val{:inf}}, point::Union{CartesianPoint{T}, PlanarPoint{T}}, ::Val{:point}) where {T}
    if point in l
        v = get_line_vector(l)
        return Line{T}[Ray(T, point, point + v), Ray(T, point, point - v)]
    else
        return Line{T}[l]
    end
end

function cut(l::Line{T,<:Any,Val{:ray}}, point::Union{CartesianPoint{T}, PlanarPoint{T}}, ::Val{:point}) where {T}
    if point == l.p1
        return Line{T}[l]
    else
        if point in l
            v = get_line_vector(l)
            return Line{T}[LineSegment(T, l.p1, point), Ray(T, point, point + v)]
        else
            return Line{T}[l]
        end
    end
end

function cut(l::Line{T,<:Any,Val{:seg}}, point::Union{CartesianPoint{T}, PlanarPoint{T}}, ::Val{:point}) where {T}
    if point == l.p1 || point == l.p2
        return Line{T}[l]
    else
        if point in l
            return Line{T}[ LineSegment(T, l.p1, point), LineSegment(T, point, l.p2)]
        else
            return Line{T}[l]
        end
    end
end

function merge(l1::Line{T,<:Any,Val{:seg}}, l2::Line{T,<:Any,Val{:seg}}) where {T}
    if isapprox(dot(get_line_vector(l1), get_line_vector(l2)), T(1), atol = geom_atol_zero(T))
        if l2.p1 in l1
            LineSegment(T, l1.p1, l2.p2), true
        elseif l2.p2 in l1
            LineSegment(T, l2.p1, l1.p2), true
        else
            return l1, false
        end
    elseif isapprox(dot(get_line_vector(l1), get_line_vector(l2)), T(-1), atol = geom_atol_zero(T))
        if l2.p1 in l1
            LineSegment(T, l1.p2, l2.p2), true
        elseif l2.p1 in l1
            LineSegment(T, l2.p1, l1.p1), true
        else
            return l1, false
        end
    else
        return l1, false
    end
end

function sample(l::Line{T,<:Any,Val{:seg}}, step::AbstractFloat) where {T}
    L = norm(l.p2 - l.p1)
    lvec = get_line_vector(l)
    return [l.p1 + δ*lvec for δ in 0:step:L]
end

function sample(l::Line{T,<:Any,Val{:seg}}, Nsamps::Int) where {T}
    if Nsamps ≤ 1
        return l.p1
    else
        L = norm(l.p2 - l.p1)
        lvec = get_line_vector(l)
        return [l.p1 + δ*lvec for δ in range(0, L, length = Nsamps)]
    end
end

function get_nodes(l::Line{T,<:Any,Val{:seg}}, step::AbstractFloat) where {T}
    nodes = sample(l, step)
    if nodes[end] ≠ l.p2
        push!(nodes, l.p2)
    end
    nodes
end

perimeter(l::Line{T,<:Any,Val{:seg}}) where {T} = T(norm(l.p2 - l.p1))
