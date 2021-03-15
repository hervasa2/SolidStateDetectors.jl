function get_cut_lines(l1::LineSegment{T}, l2::LineSegment{T}) where {T}
    sol = get_intersection(Line(l1), Line(l2))
    n_sol = _get_number_of_intersections(sol)
    if n_sol == 1
        return sol[1] in l2 ? cut(l1, sol[1], Val(:point)) : [l1]
    else
        return [l1]
    end
end

function get_cut_lines(l::LineSegment{T,PlanarPoint{T}}, c::Arc{T}) where {T}
    sol1 = get_intersection(Line(l), Circle(c))
    if !isnothing(sol1)
        sol = filter(p -> _in_angular_interval_closed(get_α_at_u_v(c, p...), c.α, tol = geom_atol_zero(T)), sol1)
        return cut(l, sol, Val(:point))
    else
        return [l]
    end
end

function get_cut_lines(c::Arc{T}, l::LineSegment{T,PlanarPoint{T}}) where {T}
    sol1 = get_intersection(Line(l), Circle(c))
    if !isnothing(sol1)
        sol = filter(p -> p in l, sol1)
        return cut(c, map(p -> get_α_at_u_v(c, p...), sol), Val(:α))
    else
        return [c]
    end
end

function get_cut_lines(c1::Arc{T}, c2::Arc{T}) where {T}
    sol1 = get_intersection(Circle(c1), Circle(c2))
    n_sol = _get_number_of_intersections(sol1)
    if n_sol > 0 && !isinf(n_sol)
        sol = filter(p -> _in_angular_interval_closed(get_α_at_u_v(c2, p...), c2.α, tol = geom_atol_zero(T)), sol1)
        return cut(c1, map(p -> get_α_at_u_v(c1, p...), sol), Val(:α))
    else
        return [c1]
    end
end

#allows ConstructiveSurfaces to use machinery of ConstructiveGeometry
get_cut_surfaces(a::AbstractLinePrimitive, b::AbstractLinePrimitive) = get_cut_lines(a,b)
