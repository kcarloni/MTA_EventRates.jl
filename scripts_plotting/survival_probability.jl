
project_dir = (@__DIR__) * "/../"
include( project_dir * "scripts/setup.jl" )

f_Psurv = load_P_survival()

begin
    fig = Figure(size=(6inch, 3inch))

    DOM_depth = 1.6km

    # prop. distance
        ax = Axis( fig[1,1], 
            xlabel="zenith θ [°]", ylabel="prop. distance [km]",
            xticks=0:15:180,
            yscale=log10,
            yminorticksvisible=true, yminorticks=IntervalsBetween(9)
        )


        θs = (0:90) * 1°
        x = calc_in_earth_prop_distance.( θs, DOM_depth )
        scatterlines!( ax, θs/1°, x/1km; 
            markersize=5, color=:black 
        )
    #


    # survival prob.
    ax = Axis( fig[1,2],
        xlabel="zenith θ [°]",
        ylabel="survival probability",
        xticks=0:15:180,
    )
    for (i, E) in enumerate( (400GeV, 500GeV, 600GeV, 700GeV, 1TeV, 5TeV, 10TeV) )
        scatterlines!( ax, θs/1°, f_Psurv.( x, E, );
            markersize=5,
            color=distinct_sequential[12][i],
            label="E = $(round(E/1TeV, sigdigits=2)) TeV"
        )
    end

    Legend( fig[0,:], ax, framevisible=false, orientation=:horizontal, nbanks=2 )

    save( project_dir * "figures/survival_prob.png", fig )
    fig
end
