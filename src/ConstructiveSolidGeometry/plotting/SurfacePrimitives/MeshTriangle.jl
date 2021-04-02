struct MeshTriangle{T}
    p1::PlanarPoint{T}
    p2::PlanarPoint{T}
    p3::PlanarPoint{T}
end

function (==)(tri1::MeshTriangle, tri2::MeshTriangle) #asumes real triangle i.e. unique vertices
    v = (tri2.p1, tri2.p2, tri2.p3)
    tri1.p1 in v && tri1.p2 in v && tri1.p3 in v
end

get_spanning_vectors(tri::MeshTriangle{T}) where {T} = (PlanarVector{T}(tri.p2 - tri.p1), PlanarVector{T}(tri.p3 - tri.p1))

function get_planar_coordinates(point::PlanarPoint, tri::MeshTriangle{T})::Tuple{T,T} where {T}
    v1, v2 = get_spanning_vectors(tri)
    v = point - tri.p1
    A = hcat(v1,v2)
    A_T = transpose(A)
    x = inv(A_T*A)*(A_T*v)
    x[1], x[2] #in units of spaning vectors
end

get_vertices(tri::MeshTriangle) = (tri.p1, tri.p2, tri.p3)

function in(point::PlanarPoint, tri::MeshTriangle{T}) where {T}
    v1, v2 = get_spanning_vectors(tri)
    if cross(v1,v2) ≠ 0
        u, v = get_planar_coordinates(point, tri)
        u ≥ T(0) && v ≥ T(0) && u+v ≤ T(1)
    else
        false
    end
end

function get_circumcircle(tri::MeshTriangle{T}) where {T}
    A = hcat(transpose(hcat(tri.p1,tri.p2,tri.p3)),SVector{3,T}(1,1,1))
    b = SVector{3,T}(tri.p1.u^2+tri.p1.v^2, tri.p2.u^2+tri.p2.v^2, tri.p3.u^2+tri.p3.v^2) #faster than dot
    x = inv(A)b
    center = PlanarPoint{T}(x[1]/2, x[2]/2)
    r = sqrt(x[3] + center.u^2 + center.v^2)
    T(r), center
end

get_midpoint(tri::MeshTriangle) = (tri.p1 + tri.p2 + tri.p3)/3

function get_plot_points(tri::MeshTriangle{T}) where {T <: AbstractFloat}
    plot_points = Vector{PlanarPoint{T}}[]
    vertices = (tri.p1, tri.p2, tri.p3)
    for i in 1:length(vertices)
        push!(plot_points, Vector{PlanarPoint{T}}([vertices[i], vertices[i%length(vertices)+1]]))
    end
    plot_points
end

@recipe function f(tri::MeshTriangle)
    seriescolor --> :orange
    linewidth --> 2
    for points in get_plot_points(tri)
        @series begin
            label := ""
            points
        end
    end
end

@recipe function f(tris::Array{<:MeshTriangle})
    seriescolor --> :orange
    linewidth --> 2
    for tri in tris
        @series begin
            label := ""
            tri
        end
    end
end

function in_circumcircle(point::PlanarPoint, tri::MeshTriangle{T}) where {T} #and not a vertex
    if point in (tri.p1, tri.p2, tri.p3)
        return false
    else
        r, center = get_circumcircle(tri)
        #println(geom_round(norm(point - center)) , " ",  geom_round(r))
        return r - norm(point - center) > geom_atol_zero(T)
    end
end

function is_adjacent_to_segment(tri::MeshTriangle, l::Line)
    v = (tri.p1, tri.p2, tri.p3)
    l.p1 in v && l.p2 in v
end

function get_3rd_vertex(tri::MeshTriangle, p1::PlanarPoint, p2::PlanarPoint)
    if tri.p1 ≠ p1 && tri.p1 ≠ p2
        return tri.p1
    elseif tri.p2 ≠ p1 && tri.p2 ≠ p2
        return tri.p2
    elseif tri.p3 ≠ p1 && tri.p3 ≠ p2
        return tri.p3
    end
end

function initialize_mesh(nodes::Array{<:PlanarPoint{T}}) where {T}
    u = map(p -> p.u, nodes)
    v = map(p -> p.v, nodes)
    uMin, uMax = minimum(u) - 0.1, maximum(u) + 0.1
    vMin, vMax = minimum(v) - 0.1, maximum(v) + 0.1
    p1 = PlanarPoint{T}(uMin, vMin)
    p2 = PlanarPoint{T}(uMax, vMin)
    p3 = PlanarPoint{T}(uMax, vMax)
    p4 = PlanarPoint{T}(uMin, vMax)
    return [p1, p2, p3, p4], [MeshTriangle{T}(p1, p2, p3), MeshTriangle{T}(p3, p4, p1)]
end

function _devide_triangle_on_edge(p::PlanarPoint, tri::MeshTriangle{T}) where {T}
    if p in Line(T, tri.p1, tri.p2)
        return [MeshTriangle{T}(p, tri.p3, tri.p1), MeshTriangle{T}(p, tri.p3, tri.p2)]
    elseif p in Line(T, tri.p2, tri.p3)
        return [MeshTriangle{T}(p, tri.p1, tri.p2), MeshTriangle{T}(p, tri.p1, tri.p3)]
    elseif p in Line(T, tri.p3, tri.p1)
        return [MeshTriangle{T}(p, tri.p2, tri.p3), MeshTriangle{T}(p, tri.p2, tri.p1)]
    else
        @warn "Point on edge of triangle contradiction"
    end
end

function divide_triangles(p::PlanarPoint, triangles::Array{<:MeshTriangle{T}}) where {T} #assumes point in triangle
    ntri = length(triangles)
    if ntri == 1
        tri = triangles[1]
        return [MeshTriangle{T}(p, tri.p1, tri.p2), MeshTriangle{T}(p, tri.p2, tri.p3), MeshTriangle{T}(p, tri.p3, tri.p1)]
    elseif ntri == 2 #p in edge between triangles
        return append!(_devide_triangle_on_edge(p, triangles[1]), _devide_triangle_on_edge(p, triangles[2]))
    else
        @warn "Node is in $ntri triangles. This should not happen!"
    end
end


function update_mesh!(nodes::Array{<:PlanarPoint}, triangles::Array{<:MeshTriangle{T}}, newnode::PlanarPoint) where {T}
    triangles_copy = deepcopy(triangles)
    idx = findall(t -> newnode in t, triangles)
    tri_candidates = divide_triangles(newnode, triangles[idx])
    deleteat!(triangles, idx)
    newtriangles = MeshTriangle{T}[]
    while length(tri_candidates) > 0
        #println(tri_candidates[1])
        if isnothing(findfirst(n -> in_circumcircle(n, tri_candidates[1]), nodes))
            #println("is good ", length(tri_candidates))
            push!(newtriangles, tri_candidates[1])
            deleteat!(tri_candidates, 1)
        else
            #println("is bad ", length(tri_candidates))
            tri = tri_candidates[1]
            iad = findfirst(t -> is_adjacent_to_segment(t, LineSegment(T, tri.p2, tri.p3)), triangles)
            if !isnothing(iad)
                adjacent_triangle = triangles[iad]
                #println("adjacent ", adjacent_triangle)
                p1 = get_3rd_vertex(adjacent_triangle, tri.p2, tri.p3)
                deleteat!(triangles, iad)
                deleteat!(tri_candidates, 1)
                s1 = MeshTriangle{T}(tri.p1, tri.p2, p1)
                s2 = MeshTriangle{T}(tri.p1, tri.p3, p1)
                prepend!(tri_candidates, [s1, s2]) #swapping diagonal
            else
                @warn "Node could not be inserted"
                break
            end
        end
    end
    push!(nodes, newnode), append!(triangles, newtriangles)
end

is_inside_triangle(tri::MeshTriangle, cps::AbstractConstructivePlanarSurface) = get_midpoint(tri) in cps

function trimesh(cps::AbstractConstructivePlanarSurface, n_arc::Real)
    nodes_object = get_nodes(cps, n_arc)
    nodes, triangles = initialize_mesh(nodes_object)
    for node in nodes_object
        update_mesh!(nodes, triangles, node)
    end
    return nodes_object, filter(t -> is_inside_triangle(t, cps), triangles)
end

struct TriMesh{T}
    x::Array
    y::Array
    z::Array
end

function mesh(cps::AbstractConstructivePlanarSurface{T}; n = 30) where {T <: AbstractFloat}
    _, triangles = trimesh(cps, 20)
    x = []
    y = []
    z = []
    plane = Plane(cps)
    R = RotXY(0.0000001,0.0000001)
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

@recipe function f(m::TriMesh{T}) where {T}
    seriestype := :surface
    linewidth := 0
    seriescolor --> :blue
    colorbar := false
    m.x, m.y, m.z
end
