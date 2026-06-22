using SolidStateDetectors
using Test
using Unitful

T = Float32

@testset "Test depletion estimation" begin
    sim = Simulation{T}(joinpath(@__DIR__, "test_config_files/BEGe_01.yaml"))
    timed_calculate_electric_potential!(sim, refinement_limits=0.01)
    id = SolidStateDetectors.determine_bias_voltage_contact_id(sim.detector)
    timed_calculate_weighting_potential!(sim, id, refinement_limits=0.01)
    SolidStateDetectors._adapt_weighting_potential_to_electric_potential_grid!(
        sim, id)
    U_est = timed_estimate_depletion_voltage(sim, check_for_depletion = false) # around 2600
    ΔU = 50u"V"
    # simulate over and under depletion voltage
    U₋ = U_est - ΔU
    U₊ = U_est + ΔU
    sim.detector = SolidStateDetector(sim.detector, contact_id=id, contact_potential=U₊)
    timed_calculate_electric_potential!(sim, refinement_limits=0.01, depletion_handling=true)
    undepleted = !is_depleted(sim.point_types)
    sim.detector = SolidStateDetector(sim.detector, contact_id=id, contact_potential=U₋)
    timed_calculate_electric_potential!(sim, refinement_limits=0.01, depletion_handling=true)
    depleted = is_depleted(sim.point_types)
    @test undepleted && depleted

    # Pass a searching range (with units)
    U_alt = timed_estimate_depletion_voltage(sim, U_est * 1.5, 0u"V", tolerance = 0.1u"V")
    @test abs(U_est - U_alt) < 5u"V"

    @test_throws Exception estimate_depletion_voltage(sim, -abs(U_est), abs(U_est))
    @test_throws Exception estimate_depletion_voltage(sim, -10, 0, tolerance = 20)
    @test_throws Exception estimate_depletion_voltage(sim, 0u"kg", 20u"kg")
    @test_logs (:info,) (:info,) (:warn, r".*not in the specified range.*") estimate_depletion_voltage(sim, U_est/3, 0)

    # `adapt_to_depletion_and_bias_voltage!` rescales the impurity density and swaps in a new
    # contact potential so the simulation matches a target depletion voltage `dep` and
    # bias voltage `bias`, without re-solving the field. Check the round-trip: after adapting,
    # the estimated depletion voltage should be ≈ `dep` and the bias contact should sit at `bias`.
    dep_target = 2000u"V"
    bias_target  = 2500u"V"
    imp_model_before = sim.detector.semiconductor.impurity_density_model
    adapt_to_depletion_and_bias_voltage!(sim, dep_target, bias_target, check_for_depletion = false, verbose = false)
    @test sim.detector.contacts[id].potential == SolidStateDetectors._parse_value(T, bias_target, SolidStateDetectors.internal_voltage_unit)
    dep_sim = estimate_depletion_voltage(sim, check_for_depletion = false, verbose = false)
    @test abs(dep_sim - dep_target) < 5u"V"
    @test sim.detector.semiconductor.impurity_density_model != imp_model_before

    # Re-run simulation in place and check depletion voltage matches again. This is a check that impurity_density_model and 
    # contact_potential where adapted correctly
    timed_calculate_electric_potential!(sim, refinement_limits=0.01, depletion_handling=true)
    @test abs(estimate_depletion_voltage(sim, check_for_depletion = false, verbose = false) - dep_sim) < 5u"V"

    # Finally, compare to fresh simulation which is changed manually 
    sim_fresh = Simulation{T}(joinpath(@__DIR__, "test_config_files/BEGe_01.yaml"))
    sim_fresh.detector = SolidStateDetector(sim_fresh.detector, contact_id = id, contact_potential = bias_target)
    sim_fresh.detector = SolidStateDetector(sim_fresh.detector, sim.detector.semiconductor.impurity_density_model)
    timed_calculate_electric_potential!(sim_fresh, refinement_limits=0.01, depletion_handling=true)
    @test abs(estimate_depletion_voltage(sim_fresh, check_for_depletion = false, verbose = false) - dep_sim) < 5u"V"
end
