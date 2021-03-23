struct CylindricalAnnulus{T,TR,TP,TZ} <: AbstractRotationalSurfacePrimitive{T}
    r::TR
    φ::TP
    z::TZ
    function CylindricalAnnulus( ::Type{T},
                   r::Union{T, <:AbstractInterval{T}},
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   z::T) where {T}
        new{T,typeof(r),typeof(φ),T}(r, φ, z)
    end
end

#Constructors
CylindricalAnnulus(c::Cone{T}; z = 0) where {T} = CylindricalAnnulus(T, get_r_at_z(c,z), c.φ, T(z))

function CylindricalAnnulus(t::Torus{T}; θ = 0) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θ = T(mod(θ,2π))
    if θ == T(0)
        rMin = t.r_torus + r_tubeMin
        rMax = t.r_torus + r_tubeMax
    elseif θ == T(π)
        rMin = t.r_torus - r_tubeMax
        rMax = t.r_torus - r_tubeMin
    else
        @error "CylindricalAnnulus not defined for torroidal cordinate θ ≠ 0 and θ ≠ π. Use ConeMantle"
    end
    r = rMin == 0 ? T(rMax) : T(rMin)..T(rMax)
    CylindricalAnnulus( T, r, t.φ, t.z)
end


function CylindricalAnnulus(; rMin = 0, rMax = 1, φMin = 0, φMax = 2π, z = 0)
    T = float(promote_type(typeof.((rMin, rMax, φMin, φMax, z))...))
    r = rMin == 0 ? T(rMax) : T(rMin)..T(rMax)
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    CylindricalAnnulus(T, r, φ, T(z))
end

CylindricalAnnulus(rMin, rMax, φMin, φMax, z) = CylindricalAnnulus(;rMin = rMin, rMax = rMax, φMin = φMin, φMax = φMax, z = z)

function CylindricalAnnulus(r::Real, z::Real)
    T = float(promote_type(typeof.((r, z))...))
    CylindricalAnnulus(T, T(r), nothing, T(z))
end

get_surface_vector(a::CylindricalAnnulus{T}) where {T} = CartesianVector{T}(0,0,1)

get_r_limits(a::CylindricalAnnulus{T, <:Union{T, AbstractInterval{T}}, <:Any}) where {T} =
    (_left_radial_interval(a.r), _right_radial_interval(a.r))

get_φ_limits(a::CylindricalAnnulus{T, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_φ_limits(a::CylindricalAnnulus{T, <:Any, <:AbstractInterval}) where {T} = (a.φ.left, a.φ.right, false)

in(p::PlanarPoint, a::CylindricalAnnulus{T, <:Any, Nothing}; check_on_plane = true) where {T} =  _in_planar_r(p, a.r)

in(p::PlanarPoint, a::CylindricalAnnulus{T, <:Any, <:AbstractInterval}; check_on_plane = true) where {T} =  _in_planar_r(p, a.r) && _in_planar_α(p, a.φ)

in(p::AbstractCoordinatePoint, a::CylindricalAnnulus{T, <:Any, Nothing}; check_on_plane = true) where {T} = _eq_z(p, a.z) && _in_cyl_r(p, a.r)

in(p::AbstractCoordinatePoint, a::CylindricalAnnulus{T, <:Any, <:AbstractInterval}) where {T} = _eq_z(p, a.z) && _in_φ(p, a.φ) && _in_cyl_r(p, a.r)

function sample(a::CylindricalAnnulus{T}, step::Real) where {T}
    rMin::T, rMax::T = get_r_limits(a)
    φMin::T, φMax::T, _ = get_φ_limits(a)
    samples = [
        CylindricalPoint{T}(r,φ,a.z)
        for r in rMin:step:rMax
        for φ in (r == 0 ? φMin : φMin:step/r:φMax)
    ]
end

function sample(a::CylindricalAnnulus{T}, Nsamps::NTuple{3,Int}) where {T}
    rMin::T, rMax::T = get_r_limits(a)
    φMin::T, φMax::T, _ = get_φ_limits(a)
    samples = [
        CylindricalPoint{T}(r,φ,a.z)
        for r in (Nsamps[1] ≤ 1 ? rMin : range(rMin, rMax, length = Nsamps[1]))
        for φ in (Nsamps[2] ≤ 1 ? φMin : range(φMin, φMax, length = Nsamps[2]))
    ]
end

function get_midpoint(a::CylindricalAnnulus{T}) where {T}
    φMin::T, φMax::T, _ = get_φ_limits(a)
    rMin::T, rMax::T = get_r_limits(a)
    sφMid::T, cφMid::T = sincos((φMax + φMin)/2)
    rMid = (rMax + rMin)/2
    CartesianPoint{T}(rMid*cφMid, rMid*sφMid, a.z)
end

function distance_to_surface(point::AbstractCoordinatePoint{T}, a::CylindricalAnnulus{T, <:Any, Nothing})::T where {T}
    point = CylindricalPoint(point)
    rMin::T, rMax::T = get_r_limits(a)
    _in_cyl_r(point, a.r) ? abs(point.z - a.z) : hypot(point.z - a.z, min(abs(point.r - rMin), abs(point.r - rMax)))
end

function distance_to_surface(point::AbstractCoordinatePoint{T}, a::CylindricalAnnulus{T, <:Any, <:AbstractInterval})::T where {T}
    pcy = CylindricalPoint(point)
    rMin::T, rMax::T = get_r_limits(a)
    φMin::T, φMax::T, _ = get_φ_limits(a)
    if _in_φ(pcy, a.φ)
        Δz = abs(pcy.z - a.z)
        return _in_cyl_r(pcy, a.r) ? Δz : hypot(Δz, min(abs(pcy.r - rMin), abs(pcy.r - rMax)))
    else
        φNear = Δ_φ(T(pcy.φ),φMin) ≤ Δ_φ(T(pcy.φ),φMax) ? φMin : φMax
        if rMin == rMax
            return norm(CartesianPoint(point)-CartesianPoint(CylindricalPoint{T}(rMin,φNear,a.z)))
        else
            return distance_to_line(CartesianPoint(point),
                                    LineSegment(T,CartesianPoint(CylindricalPoint{T}(rMin,φNear,a.z)),
                                                CartesianPoint(CylindricalPoint{T}(rMax,φNear,a.z))
                                                )
                                    )
        end
        #=Δφ = min(Δ_φ(T(point.φ),φMin),Δ_φ(T(point.φ),φMax))
        y, x = point.r .* sincos(Δφ)
        d = if x < rMin
            sqrt((rMin - x)^2 + y^2 +  Δz^2)
        elseif x > rMax
            sqrt((rMax - x)^2 + y^2 +  Δz^2)
        else
            hypot(y, Δz)
        end=#
    end
end

function cut(a::CylindricalAnnulus{T}, val::Real, ::Val{:r}) where {T}
    rMin::T, rMax::T = get_r_limits(a)
    φMin::T, φMax::T, _ = get_φ_limits(a)
    if rMin < val < rMax
        return CylindricalAnnulus{T}[
            CylindricalAnnulus(rMin, T(val), φMin, φMax, a.z),
            CylindricalAnnulus(T(val), rMax, φMin, φMax, a.z)
            ]
    else
        return CylindricalAnnulus{T}[a]
    end
end

set_φ_interval(a::CylindricalAnnulus{T}, φ::Union{Nothing, <:AbstractInterval{T}}) where {T} = CylindricalAnnulus(T, a.r, φ, a.z)

function merge(a1::CylindricalAnnulus{T}, a2::CylindricalAnnulus{T}) where {T}
    if a1 == a2
        return a1, true
    else
        if a1.φ == a2.φ && a1.z == a2.z
            rMin1::T, rMax1::T = get_r_limits(a1)
            rMin2::T, rMax2::T = get_r_limits(a2)
            φMin1::T, φMax1::T, _ = get_φ_limits(a1)
            if !isempty((rMin1..rMax1) ∩ (rMin2..rMax2))
                r = (rMin1..rMax1) ∪ (rMin2..rMax2)
                return CylindricalAnnulus(r.left, r.right, φMin1, φMax1, a1.z), true
            else
                return a1, false
            end
        elseif a1.r == a2.r && a1.z == a2.z
            union = union_angular_intervals(get_angular_interval(T, a1.φ), get_angular_interval(T, a2.φ))
            if !isnothing(union)
                φ = mod(union.right - union.left, T(2π)) == 0 ? nothing : union
                return CylindricalAnnulus(T, a1.r, φ, a1.z), true
            else
                return a1, false
            end
        else
            return a1, false
        end
    end
end

Plane(a::CylindricalAnnulus) = Plane(Val(:z), a.z)

get_decomposed_lines(a::CylindricalAnnulus{T, T, Nothing}) where {T} = AbstractLinePrimitive[Circle(a.r, PlanarPoint{T}(0,0))]

function get_decomposed_lines(a::CylindricalAnnulus{T, <:AbstractInterval{T}, Nothing}) where {T}
    rMin::T, rMax::T = get_r_limits(a)
    AbstractLinePrimitive[Circle(rMin, PlanarPoint{T}(0,0)), Circle(rMax, PlanarPoint{T}(0,0))]
end

function get_decomposed_lines(a::CylindricalAnnulus{T, T, <:AbstractInterval{T}}) where {T}
    φMin::T, φMax::T, _ = get_φ_limits(a)
    sφMin, cφMin = sincos(φMin)
    sφMax, cφMax = sincos(φMax)
    AbstractLinePrimitive[
                          Arc(T, a.r, PlanarPoint{T}(0,0), a.φ),
                          LineSegment(T, PlanarPoint{T}(0,0), PlanarPoint{T}(a.r*cφMin, a.r*sφMin)),
                          LineSegment(T, PlanarPoint{T}(0,0), PlanarPoint{T}(a.r*cφMax, a.r*sφMax))
                          ]
end

function get_decomposed_lines(a::CylindricalAnnulus{T, <:AbstractInterval{T}, <:AbstractInterval{T}}) where {T}
    rMin::T, rMax::T = get_r_limits(a)
    φMin::T, φMax::T, _ = get_φ_limits(a)
    sφMin, cφMin = sincos(φMin)
    sφMax, cφMax = sincos(φMax)
    AbstractLinePrimitive[
                          Arc(T, rMin, PlanarPoint{T}(0,0), a.φ),
                          LineSegment(T, PlanarPoint{T}(rMin*cφMin, rMin*sφMin), PlanarPoint{T}(rMax*cφMin, rMax*sφMin)),
                          LineSegment(T, PlanarPoint{T}(rMin*cφMax, rMin*sφMax), PlanarPoint{T}(rMax*cφMax, rMax*sφMax)),
                          Arc(T, rMax, PlanarPoint{T}(0,0), a.φ)
                          ]
end

function LineSegment(a::CylindricalAnnulus{T}) where {T}
    rMin::T, rMax::T = get_r_limits(a)
    LineSegment(T, PlanarPoint{T}(rMin,a.z), PlanarPoint{T}(rMax,a.z))
end

get_cross_section(a::CylindricalAnnulus) = LineSegment(a)
