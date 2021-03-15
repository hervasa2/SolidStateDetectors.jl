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

get_α_at_u_v(a::Arc{T}, u::Real, v::Real) where {T} = mod(atan(v - a.center.v, u - a.center.u), 2π) #u,v are planar coordinates

get_α_limits(a::Arc{T, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_α_limits(a::Arc{T, <:Any, <:AbstractInterval}) where {T} = (a.α.left, a.α.right, false)

translate(a::Arc{T}, t::PlanarVector{T}) where {T} = Arc(T, a.r, a.center + t, a.α)

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
