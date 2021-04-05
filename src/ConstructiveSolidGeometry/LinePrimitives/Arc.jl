struct Arc{T,TP,TH} <: AbstractLinePrimitive{T}
    r::T
    center::TP
    α::TH
    function Arc( ::Type{T},
                   r::T,
                   center::PlanarPoint{T},
                   α::Union{Nothing, <:AbstractInterval{T}}) where {T}
        new{T,typeof(center),typeof(α)}(r, center, α)
    end
end

function Arc(; r = 0, center = PlanarPoint(0,0), αMin = 0, αMax = 2π)
    T = float(promote_type(typeof.((r, αMin, αMax))..., eltype(center)))
    α = mod(T(αMax) - T(αMin), T(2π)) == 0 ? nothing : T(αMin)..T(αMax)
    Arc(T, T(r), PlanarPoint{T}(center), α)
end

Arc(r, center, αMin, αMax) = Arc(; r = r, center = center, αMin = αMin, αMax = αMax)

Circle(r::T, center::PlanarPoint{T}) where {T} = Arc(T, r, center, nothing)
Circle(a::Arc{T}) where {T} = Arc(T, a.r, a.center, nothing)

get_surface_vector(a::Arc{T}, point::PlanarPoint) where {T} =  normalize(PlanarVector{T}(point - a.center))

get_α_at_u_v(a::Arc{T}, u::Real, v::Real) where {T} = mod(atan(v - a.center.v, u - a.center.u), 2π) #u,v are planar coordinates

get_α_limits(a::Arc{T, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_α_limits(a::Arc{T, <:Any, <:AbstractInterval}) where {T} = (a.α.left, a.α.right, false)

translate(a::Arc{T}, t::PlanarVector{T}) where {T} = Arc(T, a.r, a.center + t, a.α)

function distance_to_line(point::PlanarPoint{T}, a::Arc{T, <:Any, Nothing}) where {T}
    p_t = point - a.center
    r = hypot(p_t.u, p_t.v)
    abs(a.r - r)
end

function distance_to_line(point::PlanarPoint{T}, a::Arc{T, <:Any, <:AbstractInterval}) where {T}
    αMin::T, αMax::T, _ = get_α_limits(a)
    p_t = point - a.center
    r = hypot(p_t.u, p_t.v)
    α = atan(p_t.v, p_t.u)
    if _in_angular_interval_closed(α, a.α)
        return abs(a.r - r)
    else
        αNear = Δ_φ(T(α),αMin) ≤ Δ_φ(T(α),αMax) ? αMin : αMax
        sαNear, cαNear = sincos(αNear)
        return norm(p_t - PlanarPoint{T}(a.r*cαNear, a.r*sαNear))
    end
end

function get_midpoint(a::Arc{T}) where {T}
    αMin::T, αMax::T, _ = get_α_limits(a)
    sαMid::T, cαMid::T = sincos((αMax + αMin)/2)
    PlanarPoint{T}(a.r*cαMid + a.center.u, a.r*sαMid + a.center.v)
end

function cut(a::Arc{T}, val::Real, ::Val{:α}) where {T}
    αMin::T, αMax::T, _ = get_α_limits(a)
    if _in_angular_interval_open(val, a.α)
        val_in = mod(val - αMin, T(2π)) + αMin
        return Arc{T}[
                        Arc(a.r, a.center, αMin, T(val_in)),
                        Arc(a.r, a.center, T(val_in), αMax)
                     ]
    else
        return Arc{T}[a]
    end
end

function merge(a1::Arc{T}, a2::Arc{T}) where {T}
    if a1 == a2
        return a1, true
    else
        if a1.r == a2.r && a1.center == a2.center
            α = union_angular_intervals(get_angular_interval(T, a1.α), get_angular_interval(T, a2.α))
            if !isnothing(α)
                return Arc(a1.r, a1.center, α.left, α.right), true
            else
                return a1, false
            end
        else
            return a1, false
        end
    end
end

function sample(a::Arc{T}, step::AbstractFloat) where {T}
    αMin::T, αMax::T, _ = get_α_limits(a)
    samples = [PlanarPoint{T}(a.r*cos(α)+a.center.u,a.r*sin(α)+a.center.v) for α in (a.r == 0 ? αMin : αMin:step/a.r:αMax)]
end

function sample(a::Arc{T}, Nsamps::Int) where {T}
    αMin::T, αMax::T, _ = get_α_limits(a)
    samples = [PlanarPoint{T}(a.r*cos(α)+a.center.u,a.r*sin(α)+a.center.v) for α in range(αMin, αMax, length = Nsamps)]
end

function get_nodes(a::Arc{T}, step::AbstractFloat) where {T}
    αMin::T, αMax::T, _ = get_α_limits(a)
    nodes = sample(a, step)
    endpoint = PlanarPoint{T}(a.r*cos(αMax)+a.center.u,a.r*sin(αMax)+a.center.v)
    if nodes[end] ≠ endpoint
        push!(nodes, endpoint)
    end
    nodes
end

function perimeter(a::Arc{T})::T where {T}
    αMin::T, αMax::T, _ = get_α_limits(a)
    abs(αMax - αMin)*a.r
end
