struct ToroidalAnnulus{T,TR,TB,TP,TT,TZ} <: AbstractSurfacePrimitive{T}
    r_torus::TR
    r_tube::TB
    φ::TP
    θ::TT
    z::TZ
    function ToroidalAnnulus( ::Type{T},
                   r_torus::T,
                   r_tube::Union{T, <:AbstractInterval{T}},
                   φ::T,
                   θ::Union{Nothing, <:AbstractInterval{T}},
                   z::T) where {T}
        new{T,T,typeof(r_tube),T,typeof(θ),T}(r_torus, r_tube, φ, θ, z)
    end
end

#Constructors
ToroidalAnnulus(t::Torus{T}; φ = 0) where {T} = ToroidalAnnulus( T, t.r_torus, t.r_tube, T(mod(φ,2π)), t.θ, t.z)

function ToroidalAnnulus(;r_torus = 1, r_tubeMin = 0, r_tubeMax = 1, φ = 0, θMin = 0, θMax = 2π, z = 0)
    T = float(promote_type(typeof.((r_torus, r_tubeMin, r_tubeMax, φ, θMin, θMax, z))...))
    r_tube = r_tubeMin == 0 ? T(r_tubeMax) : T(r_tubeMin)..T(r_tubeMax)
    θ = mod(T(θMax) - T(θMin), T(2π)) == 0 ? nothing : T(θMin)..T(θMax)
    ToroidalAnnulus( T, T(r_torus), r_tube, T(φ), θ, T(z))
end

ToroidalAnnulus(r_torus, r_tubeMin, r_tubeMax, φ, θMin, θMax, z) = ToroidalAnnulus(;r_torus = r_torus, r_tubeMin = r_tubeMin, r_tubeMax = r_tubeMax, φ = φ, θMin = θMin, θMax = θMax, z = z)

in(p::AbstractCoordinatePoint, t::ToroidalAnnulus{<:Any, <:Any, <:Any, <:Any, Nothing}; check_on_plane = true) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _eq_φ(p, t.φ)

in(p::AbstractCoordinatePoint, t::ToroidalAnnulus{<:Any, <:Any, <:Any, <:Any, <:AbstractInterval}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _eq_φ(p, t.φ) && _in_torr_θ(p, t.r_torus, t.θ, t.z)

get_r_tube_limits(t::ToroidalAnnulus{T}) where {T} =
    (_left_radial_interval(t.r_tube),_right_radial_interval(t.r_tube))

get_θ_limits(t::ToroidalAnnulus{T, <:Any, <:Any, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_θ_limits(t::ToroidalAnnulus{T, <:Any, <:Any, <:Any, <:AbstractInterval}) where {T} = (t.θ.left, t.θ.right, false)

function sample(t::ToroidalAnnulus{T}, step::Real) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    samples = [
        CylindricalPoint{T}(t.r_torus+r_tube*cos(θ),t.φ,r_tube*sin(θ)+t.z)
        for r_tube in r_tubeMin:step:r_tubeMax
        for θ in (r_tube == 0 ? θMin : θMin:step/r_tube:θMax)
    ]
end

function sample(t::ToroidalAnnulus{T}, Nsamps::NTuple{3,Int}) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    samples = [
        CylindricalPoint{T}(t.r_torus+r_tube*cos(θ),t.φ,r_tube*sin(θ)+t.z)
        for r_tube in (Nsamps[1] ≤ 1 ? r_tubeMin : range(r_tubeMin, r_tubeMax, length = Nsamps[1]))
        for θ in (Nsamps[3] ≤ 1 ? θMin : range(θMin, θMax, length = Nsamps[3]))
    ]
end

function Plane(t::ToroidalAnnulus{T}) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    sφ, cφ = sincos(t.φ)
    Plane(T,
          CartesianPoint{T}(t.r_torus*cφ, t.r_torus*sφ, t.z),
          CartesianPoint{T}((t.r_torus+r_tubeMax)*cφ, (t.r_torus+r_tubeMax)*sφ, t.z),
          CartesianPoint{T}(t.r_torus*cφ, t.r_torus*sφ, r_tubeMax + t.z),
          nothing
          )
end

function _transform_point_to_xy_plane_representation(point::AbstractCoordinatePoint{T}, t::ToroidalAnnulus{T})::CartesianPoint{T} where {T}
    pct = CartesianPoint(point)
    tri = Plane(t)
    CartesianPoint{T}((get_planar_coordinates(pct, tri).*norm.(get_spanning_vectors(tri)))...,distance_to_infinite_plane(pct, tri))
end

function distance_to_surface(point::AbstractCoordinatePoint{T}, t::ToroidalAnnulus{T})::T where {T}
    pct_transformed = _transform_point_to_xy_plane_representation(point, t)
    distance_to_surface(pct_transformed, CylindricalAnnulus(T, t.r_tube, t.θ, T(0)))
end

get_decomposed_lines(t::ToroidalAnnulus{T}) where {T} =
    AbstractLinePrimitive[translate(line, PlanarVector{T}(t.r_torus, t.z)) for line in get_decomposed_lines(CylindricalAnnulus(T, t.r_tube, t.θ, T(0)))]

#get_midpoint
#get_surface_vector
