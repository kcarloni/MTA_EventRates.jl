
project_dir = (@__DIR__) * "/../"
using Pkg; Pkg.activate(project_dir)

using Polynomials
using DelimitedFiles
using Interpolations

using PyCall
np = pyimport("numpy")


const R_earth = 6371km

"""
    load_muon_flux( month="december" )

Return a function ϕ_mu(E, zen). Note this flux is μ- only (does not include μ+).

The muon flux was calculated using MCEq, see `calc_surface_flux.py`.
"""
function load_muon_flux( month="december" )

    flux = np.load( project_dir * "/saved_calcs/muon_fluxes/flux_southpole_$(month).npz" )

    E_vals = flux.get("e_grid_0") * 1GeV
    zen_vals = (0:90)°
    f_mat = reduce( hcat, 
        [ (E_vals.^3) .* ( flux.get("$i") * 1/GeV/cm^2/s/sr) for i in 0:90 ]
    )

    itp = linear_interpolation(
        ( log10.(E_vals/1GeV), (zen_vals/1°) ),
        f_mat/(1GeV^2/s/sr/cm^2),
        extrapolation_bc=Line()
    )

    f_flux(E, zen) = E^(-3) * itp( log10(E/1GeV), zen/1° ) * GeV^2/cm^2/s/sr
end

"""
    calc_in_earth_prop_distance( zenith=0°, DOM_depth=1.6km; )

Calculate the distance propagated by a muon traveling in direction given by `zenith`, starting at the south pole surface. The distance is found using the law of cosines, by solving a quadratic.

(Note: `zenith = 0°` is downgoing, `zenith = 90°` is horizontal.)
"""
function calc_in_earth_prop_distance( zenith=0°, DOM_depth=1.6km; )

    poly = Polynomial([ 
        -1 * (R_earth^2 - (R_earth - DOM_depth)^2)/1km^2,
        2 * (R_earth - DOM_depth) * cos( zenith )/1km,
        1.0
    ])
    return maximum( roots(poly) ) * 1km
end

"""
    load_P_survival()

Return a function `P_surv(prop_dist, E)`, which calculates the probability a muon with inital energy `E` traveling `prop_dist` through ice will have energy greater than the muon mass (105.6 MeV) on arrival.

The survival probability was calculated using PROPOSAL, see `calc_and_save_P_surv.jl`
"""
function load_P_survival()
    P_surv = readdlm( project_dir * "saved_calcs/P_surv/mat_P.txt" )
    E_vals = vec( readdlm( project_dir * "saved_calcs/P_surv/vec_E.txt" ) * 1GeV )
    dist_vals = vec( readdlm( project_dir * "saved_calcs/P_surv/vec_d.txt" ) * 1km )

    itp = linear_interpolation(
        ( log10.(dist_vals/1km), log10.(E_vals/1GeV) ),
        P_surv,
        extrapolation_bc=Flat()
    )
    function f_Psurv( dist, E ) 
        itp( log10(dist/1km), log10(E/1GeV) )
    end
end