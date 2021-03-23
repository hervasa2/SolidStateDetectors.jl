function get_equalized_φ_planes(A::AbstractGeometry{T}, B::AbstractGeometry{T}) where {T}
    φcuts = T[]
    φ_planes = AbstractGeometry[]
    vols_B, neg = get_decomposed_volumes(B)
    append!(vols_B, neg)
    for vol in vols_B
        φMin::T, φMax::T, is_full_2π = get_φ_limits(vol)
        is_full_2π ? nothing : append!(φcuts, [φMin, φMax])
    end
    if length(φcuts) > 1
        vols_A, neg = get_decomposed_volumes(A)
        append!(vols_A, neg)
        for vol in vols_A
            for φ in φcuts
                _in_angular_interval_closed(φ, vol.φ) ? push!(φ_planes, get_cross_section(vol, φ)) : nothing
            end
        end
    end
    return φ_planes
end

function get_cut_surfaces(A::AbstractGeometry{T}, B::AbstractGeometry{T}) where {T}
    surfs_A, surfs_B = get_decomposed_surfaces(A), get_decomposed_surfaces(B)
    append!(surfs_B, get_equalized_φ_planes(B,A))
    cutsurfs_A = AbstractGeometry[]
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
