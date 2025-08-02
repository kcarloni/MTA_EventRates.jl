
include( PalomaPath * "scripts/setup_CairoMakie.jl" )
include( "../scripts/setup.jl" )

begin
# plot the prop. distance as a function of zenith, from 0° -> 180°
    fig = Figure( size=(3inch, 3inch) )
    ax = Axis( fig[1,1], 
        xlabel="zenith θ [°]", ylabel="prop. distance [km]",
        xticks=0:30:180,
        yscale=log10,
        yminorticksvisible=true, yminorticks=IntervalsBetween(9)
    )

    DOM_depth = 1.6km

    θs = (0:180) * 1°
    x = calc_in_earth_prop_distance.( θs, DOM_depth )
    scatterlines!( ax, θs/1°, x/1km; markersize=5 )

    # save( "mta_rates_julia/figures/prop_dist.png", fig )
    fig
end

begin
# plot the prop. distance as a function of zenith, from 0° -> 90°
    fig = Figure( size=(4inch, 4inch) )
    ax = Axis( fig[1,1], 
        xlabel="zenith θ [°]", ylabel="prop. distance [km]",
        xticks=0:10:90,
        # yscale=log10,
        # yminorticks=IntervalsBetween(9),
        yminorticksvisible=true, 
        yticks=0:50:160, yminorticks=0:10:160,
    )

    DOM_depth = 1.6km

    θs = (0:90) * 1°
    x = calc_in_earth_prop_distance.( θs, DOM_depth )
    scatterlines!( ax, θs/1°, x/1km; markersize=5 )

    save( "mta_rates_julia/figures/prop_dist.png", fig )
    fig
end

