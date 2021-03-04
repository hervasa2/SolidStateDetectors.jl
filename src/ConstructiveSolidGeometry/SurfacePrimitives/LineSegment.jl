function distance_to_line_segment(point::AbstractCoordinatePoint{T}, seg::Tuple{AbstractCoordinatePoint{T},AbstractCoordinatePoint{T}})::T where {T}
    point = CartesianPoint(point)
    seg = (CartesianPoint(seg[1]), CartesianPoint(seg[2]))
    v12 = normalize(CartesianVector{T}(seg[2] - seg[1]))
    v_point_1 = CartesianVector{T}(point - seg[1])
    proj_on_v12 = dot(v12,v_point_1)
    if proj_on_v12 ≤ T(0)
        return norm(seg[1] - point)
    else
        v_point_2 = CartesianVector{T}(point - seg[2])
        if dot(v12,v_point_2) ≥ T(0)
            return norm(seg[2] - point)
        else
            return sqrt(abs(dot(v_point_1,v_point_1) - proj_on_v12^2))
        end
    end
end

function distance_to_infinite_line_2D(point::SVector{2,T}, seg::Tuple{SVector{2,T},SVector{2,T}})::T where {T}
    v12 = normalize(seg[2] - seg[1])
    v_point_1 = point - seg[1]
    proj_on_v12 = dot(v12,v_point_1)
    return sqrt(abs(dot(v_point_1,v_point_1) - proj_on_v12^2))
end

function get_Ax_By_C_line_reprecentation(seg::Tuple{SVector{2,T},SVector{2,T}}) where {T}
    # Ax + By = C where C = 0 or C = 1, points given are the same will assume one is (0,0)
    M = transpose(hcat(seg[1],seg[2]))
    if det(M) == 0 #non trivial solution to homogeneus eq Mx = 0
        if seg[1][1] == 0
            return seg[2][1] == 0 ? (T(1), T(0), T(0)) : (-seg[2][2]/seg[2][1], T(1), T(0))
        else
            return -seg[1][2]/seg[1][1], T(1), T(0)
        end
    else
        b = SVector{2,T}(1,1) #solving Mx = b
        x = inv(M)*b
        return x[1], x[2], T(1)
    end
end

function get_intersection_infinite_lines(seg1::Tuple{SVector{2,T},SVector{2,T}}, seg2::Tuple{SVector{2,T},SVector{2,T}}) where {T}
    A1, B1, C1 = get_Ax_By_C_line_reprecentation(seg1)
    A2, B2, C2 = get_Ax_By_C_line_reprecentation(seg2)
    M = SMatrix{2,2,T}(A1, A2, B1, B2)
    if det(M) == 0
        if B1 == 0 # -> A1 ≠ 0 (or else would not be a line) -> B2 = 0 (since det(M) = A1B2 - A2B1 = 0) -> A2 ≠ 0
            return C1/A1 == C2/A2 ? SVector{2,T}(Inf,Inf) : SVector{2,T}(NaN,NaN)
        else
            return C1/B1 == C2/B2 ? SVector{2,T}(Inf,Inf) : SVector{2,T}(NaN,NaN)
        end
    else
        b = SVector{2,T}(C1,C2) #solving Mx = b
        return inv(M)*b
    end
end
