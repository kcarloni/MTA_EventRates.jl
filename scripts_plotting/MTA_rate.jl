
include( "../scripts/setup.jl" )

const kHz = u"kHz"

f_flux = load_muon_flux("june")
# f_flux = (E, zen) -> f_muminus_flux(E, zen) + f_muplus_flux(E, zen)
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

function calc_muon_rate_in_IceCube(depth=1.45km; 
    f_flux=load_muon_flux("june"),
    f_Psurv=load_P_survival(; E_cut=100GeV )
    )

    f(E, z) = f_flux(E, z) * f_Psurv( 
        calc_in_earth_prop_distance(z, depth),
        E
    )

    inte = integrate_f_dEdΩ( f )
    return uconvert( u"kHz", inte * 1km^2 )
end

function calc_muons_per_day_per_MTA(depth=1.6km;
    f_flux=load_muon_flux("june"),
    f_Psurv=load_P_survival(; E_cut=105MeV )
    )

    f(E, z) = f_flux(E, z) * f_Psurv( 
        calc_in_earth_prop_distance(z, depth),
        E
    )

    day = (60s * 60 * 24)
    return natural( 
        integrate_f_dEdΩ( f ) * 200cm^2 * 1day )
end


# =====================================

mu_rate_IC = calc_muon_rate_in_IceCube()
mu_per_day_MTA = calc_muons_per_day_per_MTA( 2.45km )

# =====================================

include( PalomaPath * "scripts/setup_CairoMakie.jl" )

begin 
fig = Figure( size=(7inch, 2.5inch) )
ax = Axis( fig[1,1],
    xlabel="depth [km]",
    ylabel="muons / day / MTA",
    yminorticksvisible=true,
    yticks=0:2:12,
    yminorticks=0:1:60,
    yminorgridvisible=true,
    xminorticksvisible=true,
    # xminorticks=1.45:0.1:2.5,
    xminorticks=1.4:0.1:2.7,
    xminorgridvisible=true,
)
ylims!( ax, 0, 13)
# xlims!( ax, 0, 2.5)
xlims!( ax, 1.35, 2.7)
vspan!( 1.45, 2.45; color=:grey80, alpha=0.3 )

colors=distinct_sequential[12][[1,4]]
markers=(:utriangle, :circle)
linestyles=((:dash, 1.5), :solid)

for (i, E_cut) in enumerate( (105MeV, 100GeV) )
    f_Psurv = load_P_survival(; E_cut )

    month = "december"
    j = 1
    # for (j, month) in enumerate( ("december", "june") )
        f_flux = load_muon_flux( month )

        depths = (1.35:0.1:2.7)*km
        scatterlines!( ax, 
            depths/1km, 
            calc_muons_per_day_per_MTA.(depths; f_flux, f_Psurv);
            color=lighten( colors[j], (0.6,0.4)[i]),
            linestyle=linestyles[i],
            marker=markers[i],
            label=("$month, E > $(E_cut)")
        )
    # end
end
axislegend( ax, position=:rt, nbanks=1, #framevisible=false,
    orientation=:vertical
)

# Legend(
#     [
#         [ LineElement(; color=c) for c in colors ],
#         [ LineElement(; )]
#     ]

# )

save( project_dir * "figures/rate_by_depth.png", fig )
fig
end
