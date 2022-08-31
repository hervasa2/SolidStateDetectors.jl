# """
#     _adapt_weighting_potential_to_electric_potential_grid!(sim::Simulation, contact_id::Int)
# 
# Interpolates the [`WeightingPotetial`](@ref) of the [`Contact`](@ref) with id `contact_id`
# onto the [`Grid`](@ref) of the [`ElectricPotential`](@ref) and updates it until it converges.
# 
# ## Arguments 
# * `sim::Simulation{T}`: [`Simulation`](@ref) for which the [`WeightingPotential`] should be adapted.
# * `contact_id::Int`: The `id` of the [`Contact`](@ref) at which the [`WeightingPotential`] should be adapted.
# """
function _adapt_weighting_potential_to_electric_potential_grid!(sim::Simulation, contact_id::Int)
    @assert !ismissing(sim.electric_potential) "Electric potential missing"
    if ismissing(sim.weighting_potentials[contact_id]) calculate_weighting_potential!(sim, contact_id, verbose = false) end
    if sim.weighting_potentials[contact_id].grid != sim.electric_potential.grid
        sim.weighting_potentials[contact_id] = sim.weighting_potentials[contact_id][sim.electric_potential.grid]
        #apply_initial_state!(sim, WeightingPotential, contact_id, sim.electric_potential.grid)
        update_till_convergence!(sim, WeightingPotential, contact_id)
    end
end

"""
    estimate_depletion_voltage( sim::Simulation{T}, contact_id::Int, potential_range::AbstractRange; kwargs... )::T

Returns the potential (in V) needed at the [`Contact`](@ref) with id `contact_id`
to fully deplete the detector in a given [`Simulation`](@ref). For this, all other
contact potentials are set to `0` and the potential at the specified contact is
increased or decreased according to the `potential_range`. The depletion voltage
is set to the potential for which a previously undepleted detector becomes depleted.

## Arguments 
* `sim::Simulation{T}`: [`Simulation`](@ref) for which the depletion voltage should be determined.
* `contact_id::Int`: The `id` of the [`Contact`](@ref) at which the potential is applied.
* `potential_range::AbstractRange`: Range of potentials to be tested. Must be strictly positive or negative.
    
## Keywords
* `verbose::Bool = true`: Activate or deactivate additional info output. Default is `true`.

## Example 
```julia
using SolidStateDetectors
sim = Simulation(SSD_examples[:InvertedCoax])
calculate_electric_potential!(sim)
estimate_depletion_voltage(sim, 2, 1600:1:2500)
```

!!! warn
    This method only works if the initial `sim` was calculated for a case in which the detector is fully depleted.

!!! note
    The accuracy of the result depends on the precision of the initial simulation.
    
See also [`is_depleted`](@ref).
"""
function estimate_depletion_voltage(sim::Simulation{T}, contact_id::Int,
            potential_range::AbstractRange = range(extrema(broadcast(c -> c.potential, sim.detector.contacts))..., length = 1001);
            tol = T(1e-2),
            verbose = true)::T where {T <: AbstractFloat}
    
    @assert is_depleted(sim.point_types) "This method only works for fully depleted simulations. Please increase the potentials in the configuration file to a greater value."
    @assert first(potential_range) * last(potential_range) >= 0 "Please search for depletion voltage either in a fully positive or a fully negative range. Currently you were looking in the range
        from $(first(potential_range)) to $(last(potential_range))"
    #potential_range = (last(potential_range) > 0 && first(potential_range) >= 0) || (step(potential_range) < 0) ? potential_range : reverse(potential_range)
    
    if verbose
        @info "Looking for the depletion voltage applied to contact $(contact_id) "*
              "in the range $(extrema(potential_range).*u"V") in steps of $(step(potential_range)*u"V")."
    end
    
    for c in sim.detector.contacts
        if c.potential == 0 && c.id != contact_id continue end
        if ismissing(sim.weighting_potentials[c.id]) || sim.weighting_potentials[c.id].grid != sim.electric_potential.grid
            @info "The weighting potential $(c.id) need to be defined on the same grid as the electric potential. Adjusting weighting potential $(c.id) now."
            _adapt_weighting_potential_to_electric_potential_grid!(sim, c.id)
        end
    end
    
    # ϕρ is the electric potential resulting only from the impurity density
    # and setting all potentials on the contacts to zero
    ϕρ = deepcopy(sim.electric_potential.data)
    for c in sim.detector.contacts
        if c.potential == 0 continue end
        ϕρ .-= c.potential * sim.weighting_potentials[c.id].data
    end
    
    inside = findall(p -> p & bulk_bit > 0 && p & update_bit > 0, sim.point_types.data)
    
    ϕmin, ϕmax = extrema((ϕρ .+ potential_range[1] * sim.weighting_potentials[contact_id].data)[inside])
    eps = 1e-3 * sign( potential_range[length(potential_range)÷2] )
    initial_depletion::Bool = ϕmax - ϕmin < maximum((abs(potential_range[1]), eps))
    depletion_voltage::T = NaN
    start_local_search::T = 0
    
    @showprogress for U in potential_range
        ϕmin, ϕmax = extrema((ϕρ .+ T(U) * sim.weighting_potentials[contact_id].data)[inside])
        depleted = ϕmax - ϕmin < abs(U)
        if (initial_depletion && !depleted) || (!initial_depletion && depleted)
            start_local_search = T(U)
            break
        end
    end
    local_range::AbstractRange = initial_depletion ? range(first(potential_range), start_local_search, step = step(potential_range)) : range(start_local_search, last(potential_range), step = step(potential_range))
    size_data = size(sim.point_types.data)
    @info local_range
    
    if size_data[2] == 1 # Fully phi symmetric detectors
        Ix = CartesianIndex(1,0,0)
        Iz = CartesianIndex(0,0,1)
        neighbours = [Ix,-Ix,Iz,-Iz]
    else
        Ix = CartesianIndex(1,0,0)
        Iy = CartesianIndex(0,1,0)
        Iz = CartesianIndex(0,0,1)
        neighbours = [Ix,-Ix,Iy,-Iy,Iz,-Iz]
    end
    
    scale = 0
    undepleted = 0
    eps = 1e-1
    length_idx = length(inside[1])
    @showprogress for (i,idx) in enumerate(inside)
        for scale_loc in local_range
            center_pot = T(scale_loc) * sim.weighting_potentials[contact_id].data[idx] + ϕρ[idx]
            min_pot = center_pot + eps
            max_pot = center_pot - eps
            for neighbour in neighbours
                n_idx = idx + neighbour
                if all([!(n_idx[j] == 0 || n_idx[j] > size_data[j]) for j in 1:length_idx]) # Check whether one of the indices is outside of the grid
                    local_pot = T(scale_loc) * sim.weighting_potentials[contact_id].data[n_idx] + ϕρ[n_idx]
                    if local_pot < min_pot
                        min_pot = local_pot
                    elseif local_pot > max_pot
                        max_pot = local_pot
                    end
                end
            end
            min_pot -= tol
            max_pot += tol
            if (!initial_depletion && !(min_pot>center_pot || max_pot<center_pot)) # Check whether depletion condition is met if initially undepleted for this specific index (corresponding to center pot). Jump to next index if that is the case.
                #@info scale_loc, idx, min_pot, center_pot, max_pot
                if abs(scale_loc)> abs(scale)
                    scale = T(scale_loc)
                    local_range = range(scale, last(potential_range), step = step(potential_range))
                end         
                break
            elseif (initial_depletion  && (min_pot>center_pot || max_pot<center_pot))
                @info scale_loc, idx, min_pot, center_pot, max_pot
                if abs(scale_loc) > abs(scale)
                    scale = T(scale_loc)
                end         
                break
            end
            if scale_loc == local_range[end] # If local scale is iterated all the way throughout the local range, this grid point is undepleted.
                undepleted+=1
            end         
        end
        if (!initial_depletion && undepleted > 0) # Break loop as soon as one grid point is undepleted
            break
        end
    end
    if (initial_depletion && undepleted > 0)
        depletion_voltage = scale
        @warn "Detector is already depleted at the start of the specified voltage range ($(first(potential_range)))! Calculation will be faster when going from undepleted to depleted."
    elseif (!initial_depletion && undepleted == 0)
        depletion_voltage = scale
    end
    
    if verbose
        if !isnan(depletion_voltage)
            @info "The depletion voltage of the detector is ($(depletion_voltage) ± $(abs(T(step(potential_range))))) V applied to contact $(contact_id)."
        else
            @warn "The depletion voltage is not in the specified range $(extrema(potential_range).*u"V")."
        end
    end
    return depletion_voltage
end