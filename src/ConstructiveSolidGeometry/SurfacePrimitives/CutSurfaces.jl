function get_cut_surfaces(c1::CylindricalAnnulus{T, T, Nothing}, c2::CylindricalAnnulus{T, T, Nothing}) where {T}
    return c1.z == c2.z ? cut(c1, c2.r, Val(:r)) : [c1]
end

function get_cut_surfaces(c1::CylindricalAnnulus{T, <:Any, Nothing}, c2::CylindricalAnnulus{T, <:Any, Nothing}) where {T}
    if c1.z == c2.z
        rMin::T, rMax::T = get_r_limits(c2)
        c1_cut_1 = cut(c1, rMin, Val(:r))
        c1_cuts = length(c1_cut_1) > 1 ? append!(c1_cut_1[begin:1], cut(c1_cut_1[2], rMax, Val(:r))) : cut(c1_cut_1[1], rMax, Val(:r))
        return c1_cuts
    else
        return [c1]
    end
end

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
    return rMin < r_cm < rMax ? cut(m, c.z, Val(:z)) : [m]
end

function get_cut_surfaces(m1::ConeMantle{T, <:Any, Nothing}, m2::ConeMantle{T, <:Any, Nothing}) where {T}
    seg1 = _get_2D_line_segment(m1)
    seg2 = _get_2D_line_segment(m2)
    sol = get_intersection_infinite_lines(seg1,seg2)
    zMin::T, zMax::T = get_z_limits(m2)
    if isnan(sol[1])
        return [m1]
    elseif isinf(sol[1])
        m1_cut_1 = cut(m1, zMin, Val(:z))
        m1_cuts = length(m1_cut_1) > 1 ? append!(m1_cut_1[begin:1], cut(m1_cut_1[2], zMax, Val(:z))) : cut(m1_cut_1[1], zMax, Val(:z))
        return m1_cuts
    else
        return zMin ≤ sol[2] ≤ zMax ? cut(m1, sol[2], Val(:z)) : [m1]
    end
end
