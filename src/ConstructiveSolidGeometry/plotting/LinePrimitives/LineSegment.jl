@recipe function f(l::LineSegment{T}; n = 30) where {T}
    seriescolor --> :orange
    linewidth --> 2
    [l.p1, l.p2]
end
