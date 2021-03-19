include("Line.jl")
#include("LineSegment.jl")
include("Arc.jl")
include("Intersections.jl")
include("CutLines.jl")

#return default if no method is available, some line primitives need a point to return a normal cartesian vector (surface vector)
get_surface_vector(s::AbstractLinePrimitive, point::AbstractPlanarPoint) = get_surface_vector(s)

merge(l1::AbstractLinePrimitive, l2::AbstractLinePrimitive) = l1, false
