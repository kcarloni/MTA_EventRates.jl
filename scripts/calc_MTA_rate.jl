
include( "setup.jl" )

f_flux = load_muon_flux("december")
f_Psurv = load_P_survival()


function integrate_f_dEdΩ( f; dz=1°, dlogE=0.1 )
    
    z_vals = range( 0°, 90°, step=dz )
    logE_vals = range( 2.6, 4, step=dlogE )

    function calc_integrand( z, logE )

        E = exp10(logE) * 1GeV
        dE = E * log(10) * dlogE
        dΩ = 2pi * sin(z) * (dz/1rad) * sr 

        f( E, z ) * dE * dΩ 
    end

    inte = 0 * f(1TeV, 0°) * 1TeV * 1sr
    for z in z_vals, logE in logE_vals
        inte += calc_integrand(z, logE)
    end

    return inte 
end

function calc_muon_rate_in_IceCube()

    f(E, z) = f_flux(E, z) * f_Psurv( 
        calc_in_earth_prop_distance(z, 1.45km),
        E
    )

    inte = integrate_f_dEdΩ( f )
    return uconvert( u"kHz", inte * 1km^2 )
end

function calc_muons_per_day_per_MTA(; depth=1.6km)

    f(E, z) = f_flux(E, z) * f_Psurv( 
        calc_in_earth_prop_distance(z, depth),
        E
    )

    day = (60s * 60 * 24)
    return natural( integrate_f_dEdΩ( f ) * 200cm^2 * 1day )
end

# ================================

mu_rate_IC = calc_muon_rate_in_IceCube()
mu_per_day_MTA = calc_muons_per_day_per_MTA()


println( "approx. muon rate in IceCube: ", round(typeof(1.0u"kHz"), mu_rate_IC, sigdigits=3) )
println( "events / day / MTA: ", round(mu_per_day_MTA, sigdigits=3) )
