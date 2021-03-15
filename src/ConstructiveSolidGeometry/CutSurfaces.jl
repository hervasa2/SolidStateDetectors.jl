function get_cut_surfaces(A::AbstractGeometry{T}, B::AbstractGeometry{T}) where {T}
    surfs_A, surfs_B = get_decomposed_surfaces(A), get_decomposed_surfaces(B)
    cutsurfs_A = AbstractPrimitive[]
    for a in surfs_A
        cutsurfs_a = [a]
        for b in surfs_B
            cutsurfs_a_b = []
            for a_i in cutsurfs_a
                append!(cutsurfs_a_b, get_cut_surfaces(a_i, b))
            end
            cutsurfs_a = cutsurfs_a_b
        end
        append!(cutsurfs_A, cutsurfs_a)
    end
    return cutsurfs_A
end
