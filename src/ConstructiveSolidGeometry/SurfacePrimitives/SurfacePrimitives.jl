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

merge(s1::AbstractSurfacePrimitive, s2::AbstractConstructiveSurface) = s1, false

merge(s1::AbstractConstructiveSurface, s2::AbstractSurfacePrimitive) = s1, false

get_plane_φ(s::AbstractPlanarSurfacePrimitive) = s.φ #For now all all AbstractPlanarSurfacePrimitives have a field φ. This might change.

function cut(s::AbstractRotationalSurfacePrimitive{T}, val::Real, ::Val{:φ}) where {T}
    if _in_angular_interval_open(val, s.φ)
        φMin::T, φMax::T, _ = get_φ_limits(s)
        val_in::T = mod(val - φMin, T(2π)) + φMin
        return [
                    set_φ_interval(s, φMin..val_in)
                    set_φ_interval(s, val_in..φMax)
               ]
    else
        return [s]
    end
end

function sample_border(s::AbstractSurfacePrimitive{T}, sampling) where {T}
    samples = [ point
    for line in get_decomposed_lines(s)
    for point in sample(line, sampling)  ]
end
