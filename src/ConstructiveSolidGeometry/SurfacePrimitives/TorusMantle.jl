struct TorusMantle{T,TR,TB,TP,TT,TZ} <: AbstractRotationalSurfacePrimitive{T}
    r_torus::TR
    r_tube::TB
    φ::TP
    θ::TT
    z::TZ
    function TorusMantle( ::Type{T},
                   r_torus::T,
                   r_tube::T,
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   θ::Union{Nothing, <:AbstractInterval{T}},
                   z::T) where {T}
        new{T,T,T,typeof(φ),typeof(θ),T}(r_torus, r_tube, φ, θ, z)
    end
end

#Constructors
TorusMantle(t::Torus{T}; r_tube = 1) where {T} = TorusMantle( T, t.r_torus, T(r_tube), t.φ, t.θ, t.z)

function TorusMantle(;r_torus = 1, r_tube = 1, φMin = 0, φMax = 2π, θMin = 0, θMax = 2π, z = 0)
    T = float(promote_type(typeof.((r_torus, r_tube, φMin, φMax, θMin, θMax, z))...))
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    θ = mod(T(θMax) - T(θMin), T(2π)) == 0 ? nothing : T(θMin)..T(θMax)
    TorusMantle( T, T(r_torus), T(r_tube), φ, θ, T(z))
end
TorusMantle(r_torus, r_tube, φMin, φMax, θMin, θMax, z) = TorusMantle(;r_torus = r_torus, r_tube = r_tube, φMin = φMin, φMax = φMax, θMin = θMin, θMax = θMax, z = z)

function get_surface_vector(t::TorusMantle{T}, point::AbstractCoordinatePoint)::CartesianVector{T} where {T}
    pcy = CylindricalPoint(point)
    r = pcy.r - t.r_torus
    sφ::T, cφ::T = sincos(pcy.φ)
    CartesianVector{T}(r*cφ, r*sφ, pcy.z - t.z)
end

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, Nothing, Nothing}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube, t.z)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, <:AbstractInterval, Nothing}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_φ(p, t.φ)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, Nothing, <:AbstractInterval}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_torr_θ(p, t.r_torus, t.θ, t.z)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, <:AbstractInterval, <:AbstractInterval}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_φ(p, t.φ) && _in_torr_θ(p, t.r_torus, t.θ, t.z)

function get_θ_at_z(t::TorusMantle{T}, z::Real) where {T}
    ratio = (z-t.z)/t.r_tube
    if  -1 ≤ ratio ≤ 1
        θ1 = mod(asin(ratio), 2π)
        θ2 = mod(π - θ1, 2π)
        return θ1 == θ2 ? (T(θ1),) : minmax(T(θ1),T(θ2))
    end
end

get_θ_at_r_z(t::TorusMantle{T}, r::Real, z::Real) where {T} = mod(atan(z - t.z, r - t.r_torus), 2π)

get_r_at_θ(t::TorusMantle{T}, θ::Real) where {T} = t.r_torus + t.r_tube*cos(θ)

get_φ_limits(t::TorusMantle{T, <:Any, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
get_φ_limits(t::TorusMantle{T, <:Any, <:Any, <:AbstractInterval, <:Any}) where {T} = (t.φ.left, t.φ.right, false)

get_θ_limits(t::TorusMantle{T, <:Any, <:Any, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_θ_limits(t::TorusMantle{T, <:Any, <:Any, <:Any, <:AbstractInterval}) where {T} = (t.θ.left, t.θ.right, false)

function sample(t::TorusMantle{T}, step::Real) where {T}
    φMin::T, φMax::T, _ = get_φ_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    samples = [
        CylindricalPoint{T}(get_r_at_θ(t,θ),φ,t.r_tube*sin(θ)+t.z)
        for φ in (t.r_tube == 0 ? φMin : φMin:step/t.r_tube:φMax)
        for θ in (t.r_tube == 0 ? θMin : θMin:step/t.r_tube:θMax)
    ]
end

function sample(t::TorusMantle{T}, Nsamps::NTuple{3,Int}) where {T}
    φMin::T, φMax::T, _ = get_φ_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    samples = [
        CylindricalPoint{T}(get_r_at_θ(t,θ),φ,t.r_tube*sin(θ)+t.z)
        for φ in (Nsamps[2] ≤ 1 ? φMin : range(φMin, φMax, length = Nsamps[2]))
        for θ in (Nsamps[3] ≤ 1 ? θMin : range(θMin, θMax, length = Nsamps[3]))
    ]
end

function get_midpoint(t::TorusMantle{T}) where {T}
    φMin::T, φMax::T, _ = get_φ_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    φMid = (φMax + φMin)/2
    θMid = (θMax + θMin)/2
    rMid = get_r_at_θ(t, θMid)
    CartesianPoint(CylindricalPoint{T}(get_r_at_θ(t,θMid),φMid,t.r_tube*sin(θMid)+t.z))
end

function distance_to_surface(point::AbstractCoordinatePoint{T}, t::TorusMantle{T, <:Any, <:Any, Nothing})::T where {T}
    pcy = CylindricalPoint(point)
    cy_a = CylindricalAnnulus(T,t.r_tube..t.r_tube, t.θ, T(0))
    return distance_to_surface(CartesianPoint{T}(pcy.r-t.r_torus,point.z-t.z,T(0)), cy_a)
end

function distance_to_surface(point::AbstractCoordinatePoint{T}, t::TorusMantle{T, <:Any, <:Any, <:AbstractInterval})::T where {T}
    pcy = CylindricalPoint(point)
    if _in_φ(point, t.φ)
        cy_a = CylindricalAnnulus(T,t.r_tube..t.r_tube, t.θ, T(0))
        return distance_to_surface(CartesianPoint{T}(pcy.r-t.r_torus,point.z-t.z,T(0)), cy_a)
    else
        φMin::T, φMax::T, _ = get_φ_limits(t)
        φNear::T = Δ_φ(T(pcy.φ),φMin) ≤ Δ_φ(T(pcy.φ),φMax) ? φMin : φMax
        t_a = ToroidalAnnulus(T, t.r_torus, t.r_tube..t.r_tube, φNear, t.θ, t.z)
        return distance_to_surface(point, t_a)
    end
end

function cut(t::TorusMantle{T}, val::Real, ::Val{:θ}) where {T}
    φMin::T, φMax::T, _ = get_φ_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    if _in_angular_interval_open(val, t.θ)
        val_in = mod(val - θMin, T(2π)) + θMin
        return TorusMantle{T}[
                                TorusMantle(t.r_torus, t.r_tube, φMin, φMax, θMin, T(val_in), t.z),
                                TorusMantle(t.r_torus, t.r_tube, φMin, φMax, T(val_in), θMax, t.z)
                             ]
    else
        return TorusMantle{T}[t]
    end
end

function merge(t1::TorusMantle{T}, t2::TorusMantle{T}) where {T}
    if t1 == t2
        return t1, true
    else
        if t1.z == t2.z && t1.r_torus == t2.r_torus && t1.r_tube == t2.r_tube && t1.φ == t2.φ
            φMin1::T, φMax1::T, _ = get_φ_limits(t1)
            θ = union_angular_intervals(get_angular_interval(T, t1.θ), get_angular_interval(T, t2.θ))
            if !isnothing(θ)
                return TorusMantle(t1.r_torus, t1.r_tube, φMin1, φMax1, θ.left, θ.right, t1.z), true
            else
                return t1, false
            end
        #elseif mergeinphi
        else
            return t1, false
        end
    end
end

Arc(t::TorusMantle{T}) where {T} = Arc(T, t.r_tube, PlanarPoint{T}(t.r_torus,t.z), t.θ)
Circle(t::TorusMantle{T}) where {T} = Arc(T, t.r_tube, PlanarPoint{T}(t.r_torus,t.z), nothing)

get_cross_section(t::TorusMantle) = Arc(t)
