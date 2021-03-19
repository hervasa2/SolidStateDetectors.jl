struct ConalPlane{T,TR,TP,TZ} <: AbstractPlanarSurfacePrimitive{T}
    r::TR #if tupple trapezoid/triangle, or else rectangle
    φ::TP
    z::TZ
    #if r is a Tuple, the first entry refers to the r-interval at the bottom, the second one to the r-interval at the top
    function ConalPlane( ::Type{T},
                   r::Union{T, <:AbstractInterval{T}, Tuple{T,T}, Tuple{I,I}},
                   φ::T,
                   z::Union{T, <:AbstractInterval{T}}) where {T, I<:AbstractInterval{T}}
        new{T,typeof(r),T,typeof(z)}(r, φ, z)
    end
end

ConalPlane(c::Cone{T}; φ = 0) where {T} = ConalPlane( T, c.r, T(mod(φ,2π)), c.z)

function ConalPlane(;rbotMin = 0, rbotMax = 1, rtopMin = 0, rtopMax = 1, φ = 0, zMin = -1/2, zMax = 1/2)
    T = float(promote_type(typeof.((rtopMin, rtopMax, rbotMin, rbotMax, φ, zMin, zMax))...))
    c = Cone(rbotMin, rbotMax, rtopMin, rtopMax, 0, 2π, zMin, zMax)
    ConalPlane( T, c.r, T(φ), c.z)
end
ConalPlane(rbotMin, rbotMax, rtopMin, rtopMax, φ, zMin, zMax) = ConalPlane(; rbotMin = rbotMin, rbotMax = rbotMax, rtopMin = rtopMin, rtopMax = rtopMax, φ = φ, zMin = zMin, zMax = zMax)

get_r_at_z(c::ConalPlane{T}, z::Real) where {T} = get_r_at_z(Cone(T, c.r, nothing, c.z), z::Real)
get_r_limits(c::ConalPlane{T}) where {T} = get_r_limits(Cone(T, c.r, nothing, c.z))
get_z_limits(c::ConalPlane{T}) where {T} = (_left_linear_interval(c.z), _right_linear_interval(c.z))

in(p::PlanarPoint, c::ConalPlane) =
    _in_planar_v(p, _extend_number_to_symmetric_interval(c.z)) && _in_planar_u(p, _extend_number_to_zero_interval(get_r_at_z(c, p.v)))

in(p::AbstractCoordinatePoint, c::ConalPlane) =
    _in_z(p, c.z) && _eq_φ(p, c.φ) && _in_cyl_r(p, get_r_at_z(c, p.z))

#function sample(c::ConalPlane{T}, step::Quantity{<:Real, Unitful.𝐋}) where {T}
function sample(c::ConalPlane{T}, step::Real) where {T}
    zMin::T, zMax::T = get_z_limits(c)
    #step = T(ustrip(uconvert(u"m", step)))
    samples = [
        CylindricalPoint{T}(r,c.φ,z)
        for z in zMin:step:zMax
        for r in _left_radial_interval(get_r_at_z(c, z)):step:_right_radial_interval(get_r_at_z(c, z))
    ]
end

function sample(c::ConalPlane{T}, Nsamps::NTuple{3,Int}) where {T}
    zMin::T, zMax::T = get_z_limits(c)
    samples = [
        CylindricalPoint{T}(r,c.φ,z)
        for z in (Nsamps[3] ≤ 1 ? zMin : range(zMin, zMax, length = Nsamps[3]))
        for r in (Nsamps[1] ≤ 1 ? _left_radial_interval(get_r_at_z(c, z)) : range(_left_radial_interval(get_r_at_z(c, z)), _right_radial_interval(get_r_at_z(c, z)), length = Nsamps[1]))
    ]
end

function get_vertices(c::ConalPlane{T}) where {T}
    rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
    zMin::T, zMax::T = get_z_limits(c)
    sφ, cφ = sincos(c.φ)
    (CartesianPoint{T}(rbotMin * cφ, rbotMin * sφ, zMin),
    CartesianPoint{T}(rbotMax * cφ, rbotMax * sφ, zMin),
    CartesianPoint{T}(rtopMax * cφ, rtopMax * sφ, zMax),
    CartesianPoint{T}(rtopMin * cφ, rtopMin * sφ, zMax))
end

function get_vertices_2D(c::ConalPlane{T}) where {T}
    rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
    zMin::T, zMax::T = get_z_limits(c)
    sφ, cφ = sincos(c.φ)
    (PlanarPoint{T}(rbotMin, zMin),
    PlanarPoint{T}(rbotMax, zMin),
    PlanarPoint{T}(rtopMax, zMax),
    PlanarPoint{T}(rtopMin, zMax))
end

function distance_to_surface(point::AbstractCoordinatePoint{T}, c::ConalPlane{T})::T where {T}
    rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
    zMin::T, zMax::T = get_z_limits(c)
    pcy = CylindricalPoint(point)
    Δφ = pcy.φ - c.φ
    d, r_on_plane = pcy.r .* sincos(Δφ)
    if point.z ≥ zMax
        if r_on_plane ≥ rtopMax
            return hypot(d, point.z-zMax, r_on_plane-rtopMax)
        elseif r_on_plane ≤ rtopMin
            return hypot(d, point.z-zMax, rtopMin - r_on_plane)
        else
            return hypot(d, point.z-zMax)
        end
    elseif point.z ≤ zMin
        if r_on_plane ≥ rbotMax
            return hypot(d, zMin-point.z, r_on_plane-rbotMax)
        elseif r_on_plane ≤ rtopMin
            return hypot(d, zMin-point.z, rbotMin - r_on_plane)
        else
            return hypot(d, zMax-point.z)
        end
    else
        r_at_z = get_r_at_z(c, point.z)
        rMin  = _left_radial_interval(r_at_z)
        rMax = _right_radial_interval(r_at_z)
        if rMin ≤ r_on_plane ≤ rMax
            return abs(d)
        else
            line = r_on_plane ≥ rMax ? Line(T, PlanarPoint{T}(rbotMax,zMin), PlanarPoint{T}(rtopMax,zMax)) : Line(T, PlanarPoint{T}(rbotMin,zMin), PlanarPoint{T}(rtopMin,zMax))
            point = PlanarPoint{T}(r_on_plane,point.z)
            return sqrt(d^2 + distance_to_line(point, line)^2)
        end
    end
end

Plane(c::ConalPlane) = Plane(Val(:φ), c.φ)
    #=rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
    v = get_vertices(c)
    if rbotMin == rbotMax
        return Plane(T, v[2:4]..., nothing)
    elseif rtopMin == rtopMax
        return Plane(T, v[1:3]..., nothing)
    else
        return Plane(T, v...)
    end=#

function get_decomposed_lines(c::ConalPlane{T}) where {T} #plane and get_decomposed_lines have to be defined with the same origin
    rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
    v = get_vertices_2D(c)
    if rbotMin == rbotMax
        return AbstractLinePrimitive[
                                         LineSegment(T,v[2],v[3]),
                                         LineSegment(T,v[3],v[4]),
                                         LineSegment(T,v[4],v[2]),
                                     ]
    elseif rtopMin == rtopMax
        return AbstractLinePrimitive[
                                         LineSegment(T,v[1],v[2]),
                                         LineSegment(T,v[2],v[3]),
                                         LineSegment(T,v[3],v[1]),
                                     ]
    else
        return AbstractLinePrimitive[
                                         LineSegment(T,v[1],v[2]),
                                         LineSegment(T,v[2],v[3]),
                                         LineSegment(T,v[3],v[4]),
                                         LineSegment(T,v[4],v[1])
                                     ]
    end
end

#get_midpoint
#get_surface_vector
