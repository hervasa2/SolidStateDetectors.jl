function get_cut_lines(A::AbstractGeometry{T}, B::AbstractGeometry{T}) where {T}
    lines_A, lines_B = get_decomposed_lines(A), get_decomposed_lines(B)
    cutlines_A = AbstractLinePrimitive[]
    for a in lines_A
        cutlines_a = [a]
        for b in lines_B
            cutlines_a_b = []
            for a_i in cutlines_a
                append!(cutlines_a_b, get_cut_lines(a_i, b))
            end
            cutlines_a = cutlines_a_b
        end
        append!(cutlines_A, cutlines_a)
    end
    return cutlines_A
end
