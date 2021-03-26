get_plot_points(l::Line; n = 30) = get_nodes(l, n)
get_plot_points(l::Arc; n = 30) = sample(l, n)

@recipe function f(l::AbstractLinePrimitive{T}; n = 30) where {T}
    seriescolor --> :orange
    linewidth --> 2
    get_plot_points(l)
end
