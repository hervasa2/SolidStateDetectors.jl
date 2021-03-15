#note that cut() is only defines for Line and Surface Primitives
cut(surf::AbstractPrimitive, val::Nothing, cutdir::Val) = [surf]

cut(surf::AbstractPrimitive, val::Tuple{}, cutdir::Val) = [surf]

cut(surf::AbstractPrimitive, val::Tuple{T}, cutdir::Val) where {T} = cut(surf, val[1], cutdir)

function cut(surf::AbstractPrimitive, val::Tuple{T,T}, cutdir::Val) where {T}
    cutvals = minmax(val...)
    cut1 = cut(surf, cutvals[1], cutdir)
    length(cut1) > 1 ? append!(cut1[begin:1], cut(cut1[2], cutvals[2], cutdir)) : cut(cut1[1], cutvals[2], cutdir)
end

get_cut_surfaces(c1::CylindricalAnnulus{T, <:Any, Nothing}, c2::CylindricalAnnulus{T, <:Any, Nothing}) where {T} = [c1]

function get_cut_surfaces(c::CylindricalAnnulus{T, <:Any, Nothing}, m::ConeMantle{T, <:Any, Nothing}) where {T}
    zMin::T, zMax::T = get_z_limits(m)
    if zMin ≤ c.z ≤ zMax
        r_cm = get_r_at_z(m, c.z)
        return cut(c, r_cm, Val(:r))
    else
        return [c]
    end
end

function get_cut_surfaces(m::ConeMantle{T, <:Any, Nothing}, c::CylindricalAnnulus{T, <:Any, Nothing}) where {T}
    rMin::T, rMax::T = get_r_limits(c)
    r_cm = get_r_at_z(m, c.z)
    return rMin ≤ r_cm ≤ rMax ? cut(m, c.z, Val(:z)) : [m]
end

function get_cut_surfaces(m1::ConeMantle{T, <:Any, Nothing}, m2::ConeMantle{T, <:Any, Nothing}) where {T}
    sol = get_intersection(Line(m1), Line(m2))
    n_sol = _get_number_of_intersections(sol)
    if n_sol == 1
        zMin::T, zMax::T = get_z_limits(m2)
        return zMin ≤ sol[1][2] ≤ zMax ? cut(m1, sol[1][2], Val(:z)) : [m1]
    else
        return [m1]
    end
end

function get_cut_surfaces(t::TorusMantle{T, <:Any, <:Any, Nothing}, c::CylindricalAnnulus{T, <:Any, Nothing}) where {T}
    θ_t1 = get_θ_at_z(t, c.z)
    if !isnothing(θ_t1)
        rMin::T, rMax::T = get_r_limits(c)
        θ_t = filter(θ -> rMin ≤ get_r_at_θ(t, θ) ≤ rMax, θ_t1)
        return cut(t, θ_t, Val(:θ))
    else
        return [t]
    end
end

function get_cut_surfaces(c::CylindricalAnnulus{T, <:Any, Nothing}, t::TorusMantle{T, <:Any, <:Any, Nothing}) where {T}
    θ_t1 = get_θ_at_z(t, c.z)
    if !isnothing(θ_t1)
        θ_t = filter(θ -> _in_angular_interval_closed(θ, t.θ, tol = geom_atol_zero(T)), θ_t1)
        return cut(c, map(θ -> get_r_at_θ(t, θ), θ_t), Val(:r))
    else
        return [c]
    end
end


function get_cut_surfaces(t::TorusMantle{T, <:Any, <:Any, Nothing}, m::ConeMantle{T, <:Any, Nothing}) where {T}
    sol1 = get_intersection(Line(m), Circle(t))
    if !isnothing(sol1)
        zMin::T, zMax::T = get_z_limits(m)
        sol = filter(p -> zMin ≤ p[2] ≤ zMax, sol1)
        return cut(t, map(p -> get_θ_at_r_z(t, p...), sol), Val(:θ))
    else
        return [t]
    end
end

function get_cut_surfaces(m::ConeMantle{T, <:Any, Nothing}, t::TorusMantle{T, <:Any, <:Any, Nothing}) where {T}
    sol1 = get_intersection(Line(m), Circle(t))
    if !isnothing(sol1)
        sol = filter(p -> _in_angular_interval_closed(get_θ_at_r_z(t, p...), t.θ, tol = geom_atol_zero(T)), sol1)
        return cut(m, map(p -> p[2], sol), Val(:z))
    else
        return [m]
    end
end

function get_cut_surfaces(t1::TorusMantle{T, <:Any, <:Any, Nothing}, t2::TorusMantle{T, <:Any, <:Any, Nothing}) where {T}
    sol1 = get_intersection(Circle(t1), Circle(t2))
    n_sol = _get_number_of_intersections(sol1)
    if n_sol > 0 && !isinf(n_sol)
        sol = filter(p -> _in_angular_interval_closed(get_θ_at_r_z(t2, p...), t2.θ, tol = geom_atol_zero(T)), sol1)
        return cut(t1, map(p -> get_θ_at_r_z(t1, p...), sol), Val(:θ))
    else
        return [t1]
    end
end
