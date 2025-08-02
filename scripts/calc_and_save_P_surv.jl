
project_dir = (@__DIR__) * "/../"
using Pkg; Pkg.activate(project_dir)

using DelimitedFiles

using PyCall
np = pyimport("numpy")
pp = pyimport("proposal")


function get_propagator()

    pp_args = (
        particle = pp.particle.MuMinusDef(),
        target   = pp.medium.Ice(),
        interpolate = true,
        cuts = pp.EnergyCutSettings( Inf, 0.05, false )
    )

    xsec = pp.crosssection.make_std_crosssection(
        pp_args.particle,
        pp_args.target,
        pp_args.cuts,
        pp_args.interpolate
    )

    collection = pp.PropagationUtilityCollection()
    collection.displacement = pp.make_displacement(xsec, true)
    collection.interaction = pp.make_interaction(xsec, true)
    collection.time = pp.make_time(xsec, pp_args.particle, true)
    utility = pp.PropagationUtility(collection=collection)

    geometry = pp.geometry.Sphere( pp.Cartesian3D(0,0,0), 1e20 )
    density_distr = pp.density_distribution.density_homogeneous(
        pp_args.target.mass_density)

    prop = pp.Propagator(
        pp_args.particle, 
        [(geometry, utility, density_distr)]
    )
    return prop 
end

function propagate_muon( E, prop_dist; prop )

    init_state = pp.particle.ParticleState()
    init_state.energy = E/1MeV 
    init_state.position = pp.Cartesian3D(0, 0, 0)
    init_state.direction = pp.Cartesian3D(0, 0, 1)

    track = prop.propagate(init_state, prop_dist/1cm)
end

function propagate_muon_and_return_energy( E, prop_dist; prop )
    track = propagate_muon( E, prop_dist; prop )

    if track.final_state().propagated_distance*1cm < prop_dist
        return 0MeV
    else
        return track.final_state().energy * 1MeV 
    end
end

"""
    function calc_P_survival( E, prop_dist; N_try=100, E_cut=m_muon )

Calculates the survival probability, i.e. the probability that the muon has an energy greater than `E_cut` after propagating `prop_dist` through ice. 

For speed, muons with insufficient energy `E` to travel `prop_dist`, even at MIP energy loss rates (=5m/GeV or 2MeV/cm), are immediately determined to have zero survival probability.
"""
function calc_P_survival( E, prop_dist; N_try=100, E_cut=m_muon )

    # MIP approximation: dE / dx = 2MeV/cm -> dx / dE = 5m/GeV
    if ( E * (5m/GeV) < prop_dist ); 
        println("a MIP would not make it")
        return 0
    end 

    prop = get_propagator()
    N_surv = 0
    for i in 1:N_try
        N_surv += propagate_muon_and_return_energy(E, prop_dist; prop) >= E_cut
    end
    return N_surv / N_try 
end

# ==================================

# Some rough ball-parks:

# - the PDOMs will be buried ~1.6 - 2.1km below the South Pole surface.
# - a 400GeV muon will have a 0% probability of traveling 1.6km of ice. 
# - to reach a 1.6km depth at a zenith angle of 82°, a muon will travel through about 12km of ice.
# - a 10TeV muon has about a 0.5% chance of traveling through 12km of ice.

dist_vals = round.( exp10.( 0.2:0.05:1.2 ), sigdigits=3 ) * 1km

# PROPOSAL has some weird glitches (?) where the survival probability goes to zero catastrophically. just avoid evaluating at those points. 
dist_vals[5] = 2.52km


E_vals = round.( exp10.(2.6:0.2:4 ), sigdigits=3) * 1GeV 
@time P_surv = [ 
    calc_P_survival.( E, d; N_try=100_000 ) for d in dist_vals, E in E_vals
]

writedlm( project_dir * "saved_calcs/P_surv/mat_P.txt", P_surv )
writedlm( project_dir * "saved_calcs/P_surv/vec_E.txt", E_vals/1GeV )
writedlm( project_dir * "saved_calcs/P_surv/vec_d.txt", dist_vals/1km )
