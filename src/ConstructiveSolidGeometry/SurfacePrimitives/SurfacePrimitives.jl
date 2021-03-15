include("Plane.jl")
include("ConalPlane.jl")
include("CylindricalAnnulus.jl")
include("ConeMantle.jl")
include("ToroidalAnnulus.jl")
include("TorusMantle.jl")
include("CutSurfaces.jl")
include("ConstructivePlanarSurface.jl")

#return default if no method is available, some surface primitives need a point to return a normal cartesian vector (surface vector)
get_surface_vector(s::AbstractSurfacePrimitive, point::AbstractCoordinatePoint) = get_surface_vector(s)

#if no specific merge method surfaces do not need to merge
merge(s1::AbstractSurfacePrimitive, s2::AbstractSurfacePrimitive) = s1, false

#allows ConstructiveSurfaces to use machinery of ConstructiveGeometry
get_decomposed_surfaces(s::AbstractSurfacePrimitive) = get_decomposed_lines(s)
