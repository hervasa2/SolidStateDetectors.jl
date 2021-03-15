@recipe function f(a::Arc{T}; n = 30) where {T}
    seriescolor --> :orange
    linewidth --> 2
    αMin::T, αMax::T, _ = get_α_limits(a)
    [PlanarPoint{T}(a.r*cos(α)+a.center.u,a.r*sin(α)+a.center.v) for α in range(αMin, αMax, length = n)]
end
