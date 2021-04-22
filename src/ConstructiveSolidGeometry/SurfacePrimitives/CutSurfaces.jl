#note that cut() is only defines for Line and Surface Primitives
cut(surf::AbstractPrimitive, val::Nothing, cutdir::Val) = [surf]

cut(surf::AbstractPrimitive, val::Tuple{}, cutdir::Val) = [surf]

cut(surf::AbstractPrimitive, val::Tuple{T}, cutdir::Val) where {T} = cut(surf, val[1], cutdir)

function cut(surf::AbstractPrimitive, val::Tuple{T,T}, cutdir::Val) where {T}
    cutvals = minmax(val...)
    cut1 = cut(surf, cutvals[1], cutdir)
    length(cut1) > 1 ? append!(cut1[begin:1], cut(cut1[2], cutvals[2], cutdir)) : cut(cut1[1], cutvals[2], cutdir)
end

function cut(surf::AbstractPrimitive, val::Tuple{T,T}, cutdir::Union{Val{:point}, Val{:φ}, Val{:θ}, Val{:α}}) where {T}
    cut1 = cut(surf, val[1], cutdir)
    length(cut1) > 1 ? append!(cut(cut1[1], val[2], cutdir), cut(cut1[2], val[2], cutdir)) : cut(cut1[1], val[2], cutdir)
end

#=
function get_cut_surfaces(c::CylindricalAnnulus{T}, m::ConeMantle{T}) where {T}
    if is_intersection_an_interval(get_angular_interval(T, c.φ), get_angular_interval(T, m.φ))
        zMin::T, zMax::T = get_z_limits(m)
        if zMin ≤ c.z ≤ zMax
            r_cm = get_r_at_z(m, c.z)
            return cut(c, r_cm, Val(:r))
        else
            return [c]
        end
    else
        return [c]
    end
end

function get_cut_surfaces(m::ConeMantle{T}, c::CylindricalAnnulus{T}) where {T}
    if is_intersection_an_interval(get_angular_interval(T, m.φ), get_angular_interval(T, c.φ))
        rMin::T, rMax::T = get_r_limits(c)
        r_cm = get_r_at_z(m, c.z)
        return rMin ≤ r_cm ≤ rMax ? cut(m, c.z, Val(:z)) : [m]
    else
        return [m]
    end
end

function get_cut_surfaces(m1::ConeMantle{T}, m2::ConeMantle{T}) where {T}
    if is_intersection_an_interval(get_angular_interval(T, m1.φ), get_angular_interval(T, m2.φ))
        sol = get_intersection(Line(m1), Line(m2))
        n_sol = _get_number_of_intersections(sol)
        if n_sol == 1
            zMin::T, zMax::T = get_z_limits(m2)
            return zMin ≤ sol[1][2] ≤ zMax ? cut(m1, sol[1][2], Val(:z)) : [m1]
        else
            return [m1]
        end
    else
        return [m1]
    end
end

function get_cut_surfaces(t::TorusMantle{T}, c::CylindricalAnnulus{T}) where {T}
    if is_intersection_an_interval(get_angular_interval(T, t.φ), get_angular_interval(T, c.φ))
        θ_t1 = get_θ_at_z(t, c.z)
        if !isnothing(θ_t1)
            rMin::T, rMax::T = get_r_limits(c)
            θ_t = filter(θ -> rMin ≤ get_r_at_θ(t, θ) ≤ rMax, θ_t1)
            return cut(t, θ_t, Val(:θ))
        else
            return [t]
        end
    else
        return [t]
    end
end

function get_cut_surfaces(c::CylindricalAnnulus{T}, t::TorusMantle{T}) where {T}
    if is_intersection_an_interval(get_angular_interval(T, c.φ), get_angular_interval(T, t.φ))
        θ_t1 = get_θ_at_z(t, c.z)
        if !isnothing(θ_t1)
            θ_t = filter(θ -> _in_angular_interval_closed(θ, t.θ, geom_atol_zero(T)), θ_t1)
            return cut(c, map(θ -> get_r_at_θ(t, θ), θ_t), Val(:r))
        else
            return [c]
        end
    else
        return [c]
    end
end


function get_cut_surfaces(t::TorusMantle{T}, m::ConeMantle{T}) where {T}
    if is_intersection_an_interval(get_angular_interval(T, t.φ), get_angular_interval(T, m.φ))
        sol1 = get_intersection(Line(m), Circle(t))
        if !isnothing(sol1)
            zMin::T, zMax::T = get_z_limits(m)
            sol = filter(p -> zMin ≤ p[2] ≤ zMax, sol1)
            return cut(t, map(p -> get_θ_at_r_z(t, p...), sol), Val(:θ))
        else
            return [t]
        end
    else
        return [t]
    end
end

function get_cut_surfaces(m::ConeMantle{T}, t::TorusMantle{T}) where {T}
    if is_intersection_an_interval(get_angular_interval(T, m.φ), get_angular_interval(T, t.φ))
        sol1 = get_intersection(Line(m), Circle(t))
        if !isnothing(sol1)
            sol = filter(p -> _in_angular_interval_closed(get_θ_at_r_z(t, p...), t.θ, geom_atol_zero(T)), sol1)
            return cut(m, map(p -> p[2], sol), Val(:z))
        else
            return [m]
        end
    else
        return [m]
    end
end

function get_cut_surfaces(t1::TorusMantle{T}, t2::TorusMantle{T}) where {T}
    if is_intersection_an_interval(get_angular_interval(T, t1.φ), get_angular_interval(T, t2.φ))
        sol1 = get_intersection(Circle(t1), Circle(t2))
        n_sol = _get_number_of_intersections(sol1)
        if n_sol > 0 && !isinf(n_sol)
            sol = filter(p -> _in_angular_interval_closed(get_θ_at_r_z(t2, p...), t2.θ, geom_atol_zero(T)), sol1)
            return cut(t1, map(p -> get_θ_at_r_z(t1, p...), sol), Val(:θ))
        else
            return [t1]
        end
    else
        return [t1]
    end
end
=#

get_cut_surfaces(c1::CylindricalAnnulus{T}, c2::CylindricalAnnulus{T}) where {T} = [c1]

function get_cut_surfaces(s1::AbstractRotationalSurfacePrimitive{T}, s2::AbstractRotationalSurfacePrimitive{T}) where {T}
    if is_intersection_an_interval(get_angular_interval(T, s1.φ), get_angular_interval(T, s2.φ))
        lines = get_cut_lines(get_cross_section(s1), get_cross_section(s2))
        return length(lines) == 1 ? [s1] : [get_rotational_surface_from_cross_section(l, s1.φ) for l in lines]
    else
        return [s1]
    end
end

get_cut_surfaces(s1::Union{
                               AbstractPlanarSurfacePrimitive{T},
                               PlanarSurfaceDifference{T, <:Any, <:Any, Val{:φ}},
                               PlanarSurfaceIntersection{T, <:Any, <:Any, Val{:φ}}
                           },
                 s2::AbstractRotationalSurfacePrimitive{T}
                 ) where {T} = [s1]

function get_cut_surfaces(s1::Union{
                                           AbstractPlanarSurfacePrimitive{T},
                                           PlanarSurfaceDifference{T, <:Any, <:Any, Val{:φ}},
                                           PlanarSurfaceIntersection{T, <:Any, <:Any, Val{:φ}}
                                    },
                          s2::Union{
                                           AbstractPlanarSurfacePrimitive{T},
                                           PlanarSurfaceDifference{T, <:Any, <:Any, Val{:φ}},
                                           PlanarSurfaceIntersection{T, <:Any, <:Any, Val{:φ}}
                                    }
                          ) where {T}
    if s1 == s2
        return [s1]
    elseif get_plane_φ(s1) == get_plane_φ(s2)
        intersection = PlanarSurfaceIntersection(T, s1, s2, Val(:φ))
        difference = PlanarSurfaceDifference(T, s1, s2, Val(:φ))
        if length(intersection.lines) > 0 && length(difference.lines) > 0
            return [attempt_reduce_to_primitive(difference), attempt_reduce_to_primitive(intersection)]
        else
            return [s1]
        end
    else
        return [s1]
    end
end

#only check if in φ range, requires merge to reduce number of surfaces
get_cut_surfaces(s1::AbstractRotationalSurfacePrimitive{T},
                 s2::Union{
                                AbstractPlanarSurfacePrimitive{T},
                                PlanarSurfaceDifference{T, <:Any, <:Any, Val{:φ}},
                                PlanarSurfaceIntersection{T, <:Any, <:Any, Val{:φ}}
                          }
                ) where {T} = cut(s1, get_plane_φ(s2), Val(:φ))

#does all checks, but slow
#=function get_cut_surfaces(s1::AbstractRotationalSurfacePrimitive{T}, s2::AbstractPlanarSurfacePrimitive{T}) where {T}
    φ2 = get_plane_φ(s2)
    if _in_angular_interval_open(φ2, get_angular_interval(T, s1.φ))
        l1 = get_cross_section(s1)
        lines2 = get_decomposed_lines(s2)
        cuts = 0
        for l2 in lines2
            l1cuts = get_cut_lines(l1, l2)
            cuts = cuts + length(l1cuts)
        end
        return cuts == length(lines2) ? [s1] : cut(s1, φ2, Val(:φ))
    else
        return [s1]
    end
end=#
