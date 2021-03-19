include("Line.jl")
#include("LineSegment.jl")
include("Arc.jl")

@recipe function f(l::AbstractLinePrimitive{T}; n = 30) where {T}
    seriescolor --> :orange
    linewidth --> 2
    get_plot_points(l)
end
