_get_number_of_intersections(sol::Tuple) where {T} = length(sol)
_get_number_of_intersections(sol::T) where {T} = isinf(sol) ? sol : nothing
_get_number_of_intersections(::Nothing) = 0

function get_intersection(l1::Line{T,PlanarPoint{T}}, l2::Line{T,PlanarPoint{T}}) where {T}
    A1, B1, C1 = get_Ax_By_C_line_representation(l1)
    A2, B2, C2 = get_Ax_By_C_line_representation(l2)
    M = SMatrix{2,2,T}(A1, A2, B1, B2)
    if det(M) == 0 #parallel lines
        if B1 == 0 # -> A1 ≠ 0 (or else would not be a line) -> B2 = 0 (since det(M) = A1B2 - A2B1 = 0) -> A2 ≠ 0
            return C1/A1 == C2/A2 ? Inf : nothing
        else
            return C1/B1 == C2/B2 ? Inf : nothing
        end
    else
        b = PlanarPoint{T}(C1,C2) #solving Mx = b
        return (inv(M)*b,)
    end
end

function get_intersection(l::Line{T,PlanarPoint{T}}, c::Arc{T, PlanarPoint{T}, Nothing}) where {T}
    A, B, C = get_Ax_By_C_line_representation(l)
    if B == 0
        r = C/A
        Δ = c.r^2 - (r - c.center.u)^2
        if Δ > 0
            y = sqrt(Δ)
            return (PlanarPoint{T}(r, c.center.v - y), PlanarPoint{T}(r, c.center.v + y))
        elseif Δ == 0
            return (PlanarPoint{T}(r, c.center.v),)
        end
    else
        m_l = -A/B
        b_l = C/B
        a = m_l^2 + 1
        b = 2*m_l*(b_l-c.center.v) - 2c.center.u
        c = (b_l-c.center.v)^2 - c.r^2 + c.center.u^2
        if a == 0
            x = -c/b
            return (PlanarPoint{T}(x, m_l*x + b_l),)
        else
            Δ = b^2 - 4*a*c
            if Δ > 0
                x1 = (-b + sqrt(Δ))/(2a)
                x2 = (-b - sqrt(Δ))/(2a)
                return (PlanarPoint{T}(x1, m_l*x1 + b_l), PlanarPoint{T}(x2, m_l*x2 + b_l))
            elseif Δ == 0
                x = -b/(2a)
                return (PlanarPoint{T}(x, m_l*x + b_l),)
            end
        end
    end
end

get_intersection(c::Arc{T, PlanarPoint{T}, Nothing}, l::Line{T,PlanarPoint{T}}) where {T} = get_intersection(l,c)

function get_intersection(c1::Arc{T, PlanarPoint{T}, Nothing}, c2::Arc{T, PlanarPoint{T}, Nothing}) where {T}
    d = norm(c2.center - c1.center)
    if d > c1.r + c2.r || d < abs(c1.r - c2.r)
        return nothing
    elseif d == 0 && c1.r == c2.r
        return Inf
    else
        a = (c1.r^2 - c2.r^2 + d^2)/(2d)
        h = sqrt(c1.r^2 - a^2)
        point = c1.center + (a/d)*(c2.center - c1.center)
        if h == 0
            return (point,)
        else
            Δx = h*(c2.center.v - c1.center.v)/d
            Δy = h*(c2.center.u - c1.center.u)/d
            return (PlanarPoint{T}(point[1]+Δx, point[2]-Δy), PlanarPoint{T}(point[1]-Δx, point[2]+Δy))
        end
    end
end
