
include( PalomaPath * "scripts/setup_CairoMakie.jl" )
include( "../scripts/setup.jl" )

month = "december"
f_muminus_flux = load_muon_flux( month, "mu-" )
f_muplus_flux = load_muon_flux( month, "mu+" )

begin
# plot E^3 ϕ

    fig = Figure( size=(4inch, 4inch) )
    ax = Axis( fig[1,1], 
        xscale=log10, yscale=log10,
        limits=(1e1, 1e9, 1e-6, 1e0),
        xlabel=L"E [\textrm{GeV}]",
        ylabel=L"E^{3} \phi [\textrm{GeV}^{2} \textrm{cm}^{-2} \textrm{sr}^{-1} \textrm{s}^{-1}]",
        yticks=LogTicks(-6:0),
        xticks=LogTicks(2:8),
        yminorticksvisible=true, yminorticks=IntervalsBetween(9)
    )

    u_y = GeV^(2) / sr / cm^2 / s
    E = exp10.(1:0.1:9) * GeV

    # zenith: 0 = vertical, 90 = horizontal
    for (i, f_flux) in enumerate( (f_muminus_flux, f_muplus_flux) )
        for (j, θ_zen) in enumerate( (0°, 60°, 90°) )

            ϕ = f_flux.( E, θ_zen )
            lines!( ax, E/1GeV, (E.^3 .* ϕ)/u_y;
                # label="θ = $(θ_zen)",
                linestyle=(:solid, (:dashdot,1))[i],
                color=distinct_sequential[12][2j],
                label="μ$(("-","+")[i]) at θ = $(θ_zen)"
            )
        end
    end

    axislegend(ax, position=:lb, nbanks=3, orientation=:horizontal )
    Label( fig[1,1,Top()], "atm. muon flux: South Pole in $(month)")

    save( project_dir * "figures/flux.png", fig )
    fig
end