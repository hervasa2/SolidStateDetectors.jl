function get_decomposed_surfaces(g::AbstractConstructiveGeometry{T}) where {T}
    cutsurfs_A = get_cut_surfaces(g.a,g.b)
    cutsurfs_B = get_cut_surfaces(g.b,g.a)
    real_surfs = AbstractPrimitive[]
    for surf in cutsurfs_A
        is_real_surface(surf, g) ? push!(real_surfs, surf) : nothing
    end
    for surf in cutsurfs_B
        is_real_surface(surf, g) ? push!(real_surfs, surf) : nothing
    end
    return consolidate_surfaces(real_surfs)
end

function is_real_surface(surf::AbstractSurfacePrimitive{T}, g::AbstractGeometry{T}) where {T}
    #assumes surf is cut such that it is contained in g
    tol = 10*geom_atol_zero(T)
    pt = get_midpoint(surf)
    n = get_surface_vector(surf, pt)
    pt_pos = pt + tol*n
    pt_neg = pt - tol*n
    (pt_pos in g && !(pt_neg in g)) || (pt_neg in g && !(pt_pos in g))
end

#function is_real_surface(line::AbstractLinePrimitive{T}, g::AbstractConstructiveSurface{T}) where {T}

function consolidate_surfaces(surfs::Array{AbstractPrimitive})
    surfs_copy = deepcopy(surfs)
    if length(surfs_copy) > 1
        consolidated_surfs = AbstractPrimitive[]
        i = 1
        while i ≤ length(surfs_copy)
            surf = surfs_copy[i]
            j = i + 1
            while j ≤ length(surfs_copy)
                surf, merged = merge(surf, surfs_copy[j])
                if merged
                    deleteat!(surfs_copy, j)
                else
                    j = j + 1
                end
            end
            push!(consolidated_surfs, surf)
            i = i + 1
        end
        return consolidated_surfs
    else
        return surfs_copy
    end
end
