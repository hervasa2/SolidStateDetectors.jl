include("Plane.jl")
include("InfiniteBoundedPlane.jl")
include("ConalPlane.jl")
include("CylindricalAnnulus.jl")
include("ConeMantle.jl")
include("ToroidalAnnulus.jl")
include("TorusMantle.jl")
include("ConstructivePlanarSurface.jl")
include("CutSurfaces.jl")
include("DecomposedLines.jl")
include("CutLines.jl")
include("PromoteLineToSurface.jl")

#return default if no method is available, some surface primitives need a point to return a normal cartesian vector (surface vector)
get_surface_vector(s::AbstractSurfacePrimitive, point::AbstractCoordinatePoint) = get_surface_vector(s)

#if no specific merge method surfaces do not need to merge
merge(s1::AbstractSurfacePrimitive, s2::AbstractSurfacePrimitive) = s1, false

get_plane_φ(s::AbstractPlanarSurfacePrimitive) = s.φ #For now all all AbstractPlanarSurfacePrimitives have a field φ. This might change.
