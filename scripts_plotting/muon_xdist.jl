
include( PalomaPath * "scripts/setup_CairoMakie.jl" )

# itp = load_E_frac()

# begin
# fig = Figure()
# ax = Axis( fig[1,1] )

# E_in = 1TeV
# dist = 3kpc
# x = 0:0.01:1

# # for dist in (2.5kpc, 10kpc, 13kpc )
# for E_in in (500GeV, 1TeV, 5TeV )

#     y = itp.( log10(E_in/1GeV), log10(dist/1kpc), x )
#     lines!( ax, x, y )

# end

# fig
# end

binc, mat = load_E_frac()

size(mat)

begin
fig = Figure()
ax = Axis( fig[1,1] )

x = get_edges( binc[end] )[1:(end-1)]

stairs!( ax, x, mat[15,5,:], step=:post )
stairs!( ax, x, mat[15,11,:], step=:post )


fig
end