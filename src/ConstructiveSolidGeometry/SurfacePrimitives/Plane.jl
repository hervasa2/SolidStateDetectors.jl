struct Plane{T,TP,TA} <: AbstractPlanarSurfacePrimitive{T}
    p1::TP
    p2::TP
    p3::TP
    p4::TA
    function Plane( ::Type{T},
                   p1::CartesianPoint{T},
                   p2::CartesianPoint{T},
                   p3::CartesianPoint{T},
                   p4:: Union{Nothing, CartesianPoint{T}}) where {T}
        new{T, CartesianPoint{T}, typeof(p4)}(p1, p2, p3, p4)
    end
end

#Constructors
function Triangle(;p1 = CartesianPoint(0,0,0), p2 = CartesianPoint(1,0,0), p3 = CartesianPoint(0,1,0))
    T = float(promote_type(eltype.((p1, p2, p3))...))
    Plane(T, CartesianPoint{T}(p1), CartesianPoint{T}(p2), CartesianPoint{T}(p3), nothing)
end

Triangle(p1, p2, p3) = Triangle(;p1 = p1, p2 = p2, p3 = p3)

function Plane(::Val{:φ}, val::T) where {T}
    sφ, cφ = sincos(val)
    origin = CartesianPoint{T}(0, 0, 0)
    p2 = CartesianPoint{T}(cφ, sφ, 0)
    p3 = CartesianPoint{T}(0, 0, 1)
    Plane(T, origin, p2, p3, nothing)
end

function Plane(::Val{:z}, val::T) where {T}
    origin = CartesianPoint{T}(0, 0, val)
    p2 = CartesianPoint{T}(1, 0, val)
    p3 = CartesianPoint{T}(0, 1, val)
    Plane(T, origin, p2, p3, nothing)
end

get_vertices(tri::Plane{T, CartesianPoint{T}, Nothing}) where {T} = (tri.p1, tri.p2, tri.p3)

get_spanning_vectors(plane::Plane{T}) where {T} = (CartesianVector{T}(plane.p2 - plane.p1), CartesianVector{T}(plane.p3 - plane.p1))

function get_surface_vector(plane::Plane{T})::CartesianVector{T} where {T}
    normalize(cross(get_spanning_vectors(plane)...))
end

function get_surface_vector_nonunitary(plane::Plane{T})::CartesianVector{T} where {T}
    cross(get_spanning_vectors(plane)...)
end

function distance_to_infinite_plane(point::CartesianPoint{T}, plane::Plane{T})::T where {T}
    n = get_surface_vector(plane)
    v = point - plane.p1
    abs(dot(v,n))
end

function project_on_plane(point::CartesianPoint{T}, plane::Plane{T}) where {T}
    n = get_surface_vector(plane)
    v = point - plane.p1
    d = dot(v,n)
    return point - d * n, abs(d)
end

function on_infinite_plane(point::CartesianPoint{T}, plane::Plane{T})::Bool where {T}
    isapprox(distance_to_infinite_plane(point, plane), T(0), atol = geom_atol_zero(T))
end

function get_planar_coordinates(point::CartesianPoint{T}, tri::Plane{T, CartesianPoint{T}, Nothing})::Tuple{T,T} where {T}
    v1, v2 = get_spanning_vectors(tri)
    v = point - tri.p1
    #equation to solve Ax = v where x is new 2D representation of plane. We solve by hitting it with A^T from right
    A = hcat(v1,v2)
    A_T = transpose(A)
    x = inv(A_T*A)*(A_T*v)
    x[1], x[2] #in units of spaning vectors
end

#note that for this to be a cartesian point spanning vectors must be orthogonal
get_planar_point(point::CartesianPoint{T}, tri::Plane{T, CartesianPoint{T}, Nothing}) where {T} =
    PlanarPoint{T}((get_planar_coordinates(point, tri) .* norm.(get_spanning_vectors(tri)))...)

function get_cartesian_point(point::PlanarPoint{T}, tri::Plane{T, CartesianPoint{T}, Nothing}) where {T}
    u_vec, v_vec = point .* normalize.(get_spanning_vectors(tri))
    tri.p1 + u_vec + v_vec
end

function projection_in_triangle(point::CartesianPoint{T}, tri::Plane{T, CartesianPoint{T}, Nothing})::Bool where {T}
    u, v = get_planar_coordinates(point, tri)
    u ≥ T(0) && v ≥ T(0) && u+v ≤ T(1)
end

in(point::CartesianPoint{T}, tri::Plane{T, CartesianPoint{T}, Nothing}) where {T} = projection_in_triangle(point, tri) && on_infinite_plane(point, tri)
#Order of checks makes no difference in performance in all tested cases


function distance_to_surface(point::AbstractCoordinatePoint{T}, tri::Plane{T, CartesianPoint{T}, Nothing})::T where {T}
    point = CartesianPoint(point)
    u, v = get_planar_coordinates(point, tri)
    if geom_round(u) ≥ T(0) && geom_round(v) ≥ T(0) #++ quadrant. Origin is tri.p1
        if geom_round(u+v) ≤ T(1) #in triangle
            return distance_to_infinite_plane(point, tri)
        else # on side 2_3 of tri
            return distance_to_line(point, LineSegment(T,tri.p2,tri.p3))
        end
    elseif geom_round(u) ≤ T(0) && geom_round(v) ≤ T(0) #-- quadrant
        return norm(tri.p1 - point)
    elseif geom_round(u) > T(0) && geom_round(v) < T(0) #+- quadrant, on side 1_2 of tri
        return distance_to_line(point, LineSegment(T,tri.p1,tri.p2))
    elseif geom_round(u) < T(0) && geom_round(v) > T(0) #-+ quadrant, on side 3_1 of tri
        return distance_to_line(point, LineSegment(T,tri.p3,tri.p1))
    end
end

function sample(tri::Plane{T, CartesianPoint{T}, Nothing}, step::Real) where {T}
    v1, v2 = get_spanning_vectors(tri)
    step_u = step/norm(v1)
    step_v = step/norm(v2)
    samples = [
        CartesianPoint{T}(tri.p1 + u*v1 + v*v2)
        for u in 0:step_u:1
        for v in 0:step_v:1-u
    ]
end

function sample(tri::Plane{T, CartesianPoint{T}, Nothing}, Nsamps::NTuple{3,Int}) where {T}
    v1, v2 = get_spanning_vectors(tri)
    samples = [
        CartesianPoint{T}(tri.p1 + u*v1 + v*v2)
        for u in (Nsamps[1] ≤ 1 ? 0 : range(0, 1, length = Nsamps[1]))
        for v in (Nsamps[2] ≤ 1 ? 0 : range(0, 1-u, length = Nsamps[2]))
    ]
end

function Quadrilateral(;p1 = CartesianPoint(0,0,0), p2 = CartesianPoint(1,0,0), p3 = CartesianPoint(1,1,0), p4 = CartesianPoint(0,1,0), p4_on_plane_check = true)
    T = float(promote_type(eltype.((p1, p2, p3, p4))...))
    p1, p2, p3, p4 = CartesianPoint{T}(p1), CartesianPoint{T}(p2), CartesianPoint{T}(p3), CartesianPoint{T}(p4)
    #will return triangle if conditions are not met. Order of points matters, a continous non intersecting line needs to be drawn in p1->p2->p3->p4->p1
    if geom_round(p4) in [geom_round(p1), geom_round(p2), geom_round(p3)]
        return Plane(T, p1, p2, p3, nothing)
    elseif p4_on_plane_check
        tri = Plane(T, p1, p2, p3, nothing)
        if on_infinite_plane(p4, tri)
            if !projection_in_triangle(p4, tri)
                v = CartesianVector{T}(p4 - p1)
                v2 = CartesianVector{T}(p3 - p1)
                n = get_surface_vector_nonunitary(tri)
                if geom_round(dot(n,cross(v2,v))) ≥ 0
                    return Plane(T, p1, p2, p3, p4)
                else #intersecting segments
                    return tri
                end
            else
                return Plane(T, p1, p2, p3, p4)
            end
        else #not on plane
            return tri
        end
    else
        return Plane(T, p1, p2, p3, p4)
    end
end

Quadrilateral(p1, p2, p3, p4; p4_on_plane_check = true) = Quadrilateral(;p1 = p1, p2 = p2, p3 = p3, p4 = p4, p4_on_plane_check = p4_on_plane_check)
Plane(p1, p2, p3, p4; p4_on_plane_check = true) = Quadrilateral(;p1 = p1, p2 = p2, p3 = p3, p4 = p4, p4_on_plane_check = p4_on_plane_check)
Plane(p1, p2, p3; p4_on_plane_check = true) = Triangle(;p1 = p1, p2 = p2, p3 = p3)

get_vertices(quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}}) where {T} = (quad.p1, quad.p2, quad.p3, quad.p4)


function decompose_into_tiangles(quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}})::Tuple{Plane{T, CartesianPoint{T}, Nothing},Plane{T, CartesianPoint{T}, Nothing}} where {T}
    quad.p4 in Triangle(quad.p2,quad.p1,quad.p3) || quad.p2 in Triangle(quad.p4,quad.p1,quad.p3) ? (Triangle(quad.p1,quad.p2,quad.p4), Triangle(quad.p3,quad.p2,quad.p4)) : (Triangle(quad.p2,quad.p1,quad.p3), Triangle(quad.p4,quad.p1,quad.p3))
end


function in(point::CartesianPoint{T}, quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}})::Bool where {T}
    tri1, tri2 = decompose_into_tiangles(quad)
    point in tri1 || point in tri2
end

function distance_to_surface(point::AbstractCoordinatePoint{T}, quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}})::T where {T}
    tri1, tri2 = decompose_into_tiangles(quad)
    min(distance_to_surface(point, tri1), distance_to_surface(point, tri2))
end

function sample(quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}}, sampling::Union{Real, NTuple{3,Int}}) where {T}
    tri1, tri2 = decompose_into_tiangles(quad)
    samples = sample(tri1, sampling)
    append!(samples, sample(tri2, sampling))
end
