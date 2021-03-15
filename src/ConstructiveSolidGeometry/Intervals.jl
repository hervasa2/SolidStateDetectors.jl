@inline _left_linear_interval(z::Real) = -z
@inline _left_linear_interval(z::AbstractInterval) = z.left

@inline _right_linear_interval(z::Real) = z
@inline _right_linear_interval(z::AbstractInterval) = z.right

@inline _width_linear_interval(z::Real) = 2z
@inline _width_linear_interval(z::AbstractInterval) = width(z)


@inline _left_radial_interval(r::T) where {T <: Real} = T(0)
@inline _left_radial_interval(r::AbstractInterval) = r.left

@inline _right_radial_interval(r::Real) = r
@inline _right_radial_interval(r::AbstractInterval) = r.right

@inline _extend_number_to_zero_interval(r::T) where {T <: Real} = T(0)..r
@inline _extend_number_to_zero_interval(r::AbstractInterval) = r

@inline _in_angular_interval_closed(α::T, α_int::Nothing; tol = 0) where {T<:Real} = true
@inline _in_angular_interval_closed(α::Real, α_int::AbstractInterval{T}; tol = 0) where {T} = -tol ≤ mod(α - (α_int.left-tol), T(2π)) ≤ (α_int.right+tol) - (α_int.left-tol)
@inline _in_angular_interval_open(α::T, α_int::Nothing) where {T<:Real} = 0 < mod(α, T(2π)) < 2π
@inline _in_angular_interval_open(α::Real, α_int::AbstractInterval{T}) where {T} = 0 < mod(α - α_int.left, T(2π)) < α_int.right - α_int.left

_in_angular_interval_union(val::T, α::AbstractInterval, β::AbstractInterval) where {T<:AbstractFloat} =
    _in_angular_interval_closed(val, α) || _in_angular_interval_closed(val, β)

function _is_edge_of_angular_interval_union(val::T, α::AbstractInterval, β::AbstractInterval) where {T<:AbstractFloat}
    tol = 10*geom_atol_zero(T)
    (_in_angular_interval_union(val+tol, α, β) && !_in_angular_interval_union(val-tol, α, β)) ||
    (_in_angular_interval_union(val-tol, α, β) && !_in_angular_interval_union(val+tol, α, β))
 end

function union_angular_intervals(α::AbstractInterval{T}, β::AbstractInterval{T}) where {T} #if no intersection will return nothing
    if α.left < 0 || β.left < 0
        edges = filter(p -> _is_edge_of_angular_interval_union(p, α, β), (α.left, α.right, β.left, β.right))
        if length(edges) == 2
            θMin, θMax = minmax(mod.(edges, T(2π))...)
            return _in_angular_interval_union(θMin+tol, α, β) ? θMin..θMax : (θMax-T(2π))..θMin
        end
    else
        return isempty(α ∩ β) ? nothing : α ∪ β
    end
end
