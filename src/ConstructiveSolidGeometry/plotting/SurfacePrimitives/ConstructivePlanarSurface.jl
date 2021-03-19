function get_plot_points(cps::AbstractConstructivePlanarSurface{T}; n = 30) where {T <: AbstractFloat}
    plot_points = Vector{CartesianPoint{T}}[]
    plane = Plane(cps)

    for line in cps.lines
        push!(plot_points, get_cartesian_point.(get_plot_points(line, n = n), (plane,)))
    end
    plot_points
end

@recipe function f(cps::AbstractConstructivePlanarSurface{T}; n = 30) where {T}
    seriescolor --> :orange
    linewidth --> 2
    for points in get_plot_points(cps, n = n)
        @series begin
            label := ""
            points
        end
    end
end
