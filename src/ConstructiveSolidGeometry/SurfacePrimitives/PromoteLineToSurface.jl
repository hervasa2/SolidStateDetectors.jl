function get_rotational_surface_from_cross_section(l::Line{T, <:Any, Val{:seg}}, φ::Union{Nothing, <:AbstractInterval{T}}) where {T}
    if l.p1.v == l.p2.v
        rMin, rMax = minmax(l.p1.u, l.p2.u)
        r = rMin == 0 ? T(rMax) : T(rMin)..T(rMax)
        return CylindricalAnnulus(T, r, φ, l.p1.v)
    else
        ptop, pbot = l.p1.v > l.p2.v ? (l.p1, l.p2) : (l.p2, l.p1)
        r = pbot.u == ptop.u ? pbot.u : (pbot.u, ptop.u)
        z = ptop.v == -pbot.v ? ptop.v : pbot.v..ptop.v
        return ConeMantle(T, r, φ, z)
    end
end

get_rotational_surface_from_cross_section(a::Arc{T}, φ::Union{Nothing, <:AbstractInterval{T}}) where {T} =
    TorusMantle(T, a.center.u, a.r, φ, a.α, a.center.v)
