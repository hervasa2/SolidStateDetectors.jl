function get_plot_points(cps::AbstractConstructivePlanarSurface{T}; n = 30) where {T <: AbstractFloat}
    plot_points = Vector{CartesianPoint{T}}[]
    plane = Plane(cps)

    for line in cps.lines
        push!(plot_points, get_cartesian_point.(get_plot_points(line, n = n), (plane,)))
    end
    plot_points
end

function mesh(cps::Union{PlanarSurfaceDifference{T, <:Any, <:Any, Val{:φ}}, PlanarSurfaceIntersection{T, <:Any, <:Any, Val{:φ}}}; n = 30) where {T}
    _, triangles = trimesh(cps, n)
    x = []
    y = []
    z = []
    plane = Plane(cps)
    R = RotXY(0.0000001,0.0000001) #vertical meshes are not supported by plots
    for tri in triangles
        p1 = R*get_cartesian_point(tri.p1, plane)
        p2 = R*get_cartesian_point(tri.p2, plane)
        p3 = R*get_cartesian_point(tri.p3, plane)
        push!(x, [p1.x, p2.x, p3.x])
        push!(y, [p1.y, p2.y, p3.y])
        push!(z, [p1.z, p2.z, p3.z])
    end
    TriMesh{T}(x, y, z)
end
