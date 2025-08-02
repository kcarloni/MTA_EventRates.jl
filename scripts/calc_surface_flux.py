from MCEq.core import MCEqRun
import crflux.models as crf
import numpy as np

month = 'June'

# Initialize MCEq
mceq = MCEqRun(
    # High-energy hadronic interaction model
    interaction_model = 'SIBYLL23C',
    # Atmospheric density model
    density_model = ('MSIS00_IC', ('SouthPole', month)),
    # CR flux at top of atmosphere
    primary_model = (crf.HillasGaisser2012, 'H3a'),
    # zenith angle in degrees
    theta_deg = 0.
)

altitude = 2835 * 1e2 # South Pole altitude in cm
X = mceq.density_model.h2X([altitude])

# Loop over theta from 0 to 90 deg, computing flux at each
fluxes = {}
for theta in range(91):
    print(f'Computing flux for theta: {theta} deg')

    mceq.set_theta_deg(theta)
    mceq.solve(int_grid=X)

    fluxes[f'e_grid_{theta}'] = mceq.e_grid
    fluxes[str(theta)] = mceq.get_solution('mu-')

np.savez(f'data/surface_fluxes/flux_southpole_{month.lower()}.npz', **fluxes)
