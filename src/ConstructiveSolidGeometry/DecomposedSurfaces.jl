function get_decomposed_surfaces(g::AbstractConstructiveGeometry{T}) where {T}
    cutsurfs_A = get_cut_surfaces(g.a,g.b)
    cutsurfs_B = get_cut_surfaces(g.b,g.a)
    real_surfs = AbstractSurfacePrimitive[]
    for surf in cutsurfs_A
        is_real_surface(surf, g) ? push!(real_surfs, surf) : nothing
    end
    for surf in cutsurfs_B
        is_real_surface(surf, g) ? push!(real_surfs, surf) : nothing
    end
    return real_surfs
end

function is_real_surface(surf::AbstractSurfacePrimitive{T}, g::AbstractGeometry{T}) where {T}
    tol = 10*geom_atol_zero(T)
    pt = get_midpoint(surf)
    pcy = CylindricalPoint(pt)
    n = get_surface_vector(surf, pcy.φ)
    pt_pos = pt + tol*n
    pt_neg = pt - tol*n
    (pt_pos in g && !(pt_neg in g)) || (pt_neg in g && !(pt_pos in g))
end
