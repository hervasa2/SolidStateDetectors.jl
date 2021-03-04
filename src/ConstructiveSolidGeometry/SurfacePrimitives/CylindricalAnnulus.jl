struct CylindricalAnnulus{T,TR,TP,TZ} <: AbstractSurfacePrimitive{T}
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
    CylindricalAnnulus( T, r, t.φ, T(0))
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

in(p::AbstractCoordinatePoint, a::CylindricalAnnulus{T, <:Any, Nothing}) where {T} = _eq_z(p, a.z) && _in_cyl_r(p, a.r)

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
            return distance_to_line_segment(point, (CylindricalPoint{T}(rMin,φNear,a.z),CylindricalPoint{T}(rMax,φNear,a.z)))
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
        return [CylindricalAnnulus(rMin, T(val), φMin, φMax, a.z), CylindricalAnnulus(T(val), rMax, φMin, φMax, a.z)]
    else
        return [a]
    end
end
