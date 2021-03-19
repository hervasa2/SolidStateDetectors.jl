function get_decomposed_lines(g::AbstractConstructiveSurface)
    cutlines_A = get_cut_lines(g.a,g.b)
    cutlines_B = get_cut_lines(g.b,g.a)
    real_lines = AbstractLinePrimitive[]
    for line in cutlines_A
        is_real_surface(line, g) ? push!(real_lines, line) : nothing
    end
    for line in cutlines_B
        is_real_surface(line, g) ? push!(real_lines, line) : nothing
    end
    return consolidate_surfaces(real_lines)
end
