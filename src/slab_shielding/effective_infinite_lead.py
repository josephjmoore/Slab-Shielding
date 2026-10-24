import openmc
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from params import (cross_sections, S, E, MeV_to_eV, lead_density,
                    mu, air_energy, air_mu_en_rho, MeV_to_J, g_to_kg, s_to_h,
                    e_dep, gp_buildup, Gy_to_uGy)

#Create Shield Material
shield = openmc.Material(name='lead shield')
shield.add_element('Pb', 1.0)
shield.set_density('g/cm3', lead_density)
#Depths of the detector in the lead to sweep in cm
shield_thickness = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] #cm
#Lead beyond the detector in cm (~12 mfp, effectively infinite)
lead_beyond_detector = 10.0 #cm

#Conversion from MeV/g/s to Gy/h
MeV_per_g_s_to_Gy_per_h = MeV_to_J * g_to_kg * s_to_h

def build_infinite_geometry(shield_thickness, r_end=10.0):
    """Build the infinite-medium geometry for one lead thickness (cm).

    Lead fills everything from the source out to r_end (cm), with a thin
    lead-filled detector shell centred at a depth of shield_thickness.
    Returns the geometry and the detector cell for tallies.
    """
    #Thin detector shell centred on the shield thickness
    r_det_in = shield_thickness - 0.05 #cm
    r_det_out = shield_thickness + 0.05 #cm

    r_detector = shield_thickness

    V_det = (4.0 / 3.0) * np.pi * (r_det_out**3 - r_det_in**3)
    flux_no_shield = S / (4 * np.pi * r_detector**2)

    #Surfaces: one sphere per radius, shared by the cells either side of it
    detector_inner_surface = openmc.Sphere(r=r_det_in)
    detector_outer_surface = openmc.Sphere(r=r_det_out)
    #Only the outermost surface is a vacuum boundary (kills escaping particles)
    outer_boundary = openmc.Sphere(r=r_end, boundary_type='vacuum')

    shield_before_cell = openmc.Cell(name='lead shield before', region=-detector_inner_surface)
    shield_before_cell.fill = shield

    detector_cell = openmc.Cell(name='detector', region=+detector_inner_surface & -detector_outer_surface)
    detector_cell.fill = shield

    shield_after_cell = openmc.Cell(name='lead shield after', region=+detector_outer_surface & -outer_boundary)
    shield_after_cell.fill = shield

    cells = [detector_cell, shield_before_cell, shield_after_cell]

    root = openmc.Universe(cells=cells)
    geometry = openmc.Geometry(root)
    return geometry, detector_cell, V_det, flux_no_shield

#Caesium-137 source
source = openmc.IndependentSource()
source.space = openmc.stats.Point((0.0, 0.0, 0.0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([E * MeV_to_eV], [1.0]) #MeV to eV
source.particle = 'photon'

#settings
settings = openmc.Settings()
settings.source = source
settings.run_mode = 'fixed source'
settings.photon_transport = True
settings.batches = 10
settings.particles = 100000

materials = openmc.Materials([shield])
openmc.config['cross_sections'] = cross_sections


def run_infinite(shield_thickness, lead_beyond):
    """Run the infinite-medium model for one depth and return the results."""
    r_end = shield_thickness + lead_beyond
    geometry, detector_cell, V_det, flux_no_shield = build_infinite_geometry(shield_thickness, r_end)

    #Tallies
    cell_filter = openmc.CellFilter(detector_cell)
    tally = openmc.Tally(name='detector tally')
    tally.filters.append(cell_filter)
    tally.scores.append('flux')

    uncollided_tally = openmc.Tally(name='uncollided tally')
    energy_filter = openmc.EnergyFilter([E * MeV_to_eV - 100, E * MeV_to_eV + 100]) #MeV to eV
    uncollided_tally.filters.append(energy_filter)
    uncollided_tally.filters.append(cell_filter)
    uncollided_tally.scores.append('flux')

    air_kerma_tally = openmc.Tally(name='air kerma tally')
    energy_function_filter = openmc.EnergyFunctionFilter(air_energy * MeV_to_eV, air_energy * air_mu_en_rho) #MeV to eV for first param as openmc compares in ev, second param in MeV*cm^2/g
    energy_function_filter.interpolation = 'log-log'
    air_kerma_tally.filters.append(energy_function_filter)
    air_kerma_tally.filters.append(cell_filter)
    air_kerma_tally.scores.append('flux')

    tallies = openmc.Tallies([tally, uncollided_tally, air_kerma_tally])

    model = openmc.Model(geometry=geometry, materials=materials, settings=settings, tallies=tallies)

    run_dir = Path(__file__).resolve().parent / 'runs' / 'infinite' / f't_{shield_thickness}cm_beyond_{lead_beyond}cm'
    statepoint_path = model.run(cwd=run_dir)
    with openmc.StatePoint(statepoint_path) as sp:
        uncollided_result = sp.get_tally(name='uncollided tally')
        uncollided_tally_mean = uncollided_result.mean.item()
        uncollided_tally_err = uncollided_result.std_dev.item()

        dose_result = sp.get_tally(name='air kerma tally')
        dose_tally_mean = dose_result.mean.item()
        dose_tally_err = dose_result.std_dev.item()

    #Uncollided flux
    uncollided_flux_openmc = uncollided_tally_mean / V_det * S

    K_dot_openmc = dose_tally_mean / V_det * S * MeV_per_g_s_to_Gy_per_h
    K_dot_openmc_err = dose_tally_err / V_det * S * MeV_per_g_s_to_Gy_per_h

    K_dot_uncollided_openmc = uncollided_flux_openmc * e_dep * MeV_per_g_s_to_Gy_per_h

    #Monte Carlo buildup factor
    B_openmc = K_dot_openmc / K_dot_uncollided_openmc
    #Error on B, two air kerma rates in quadrature
    B_openmc_err = B_openmc * np.sqrt((dose_tally_err / dose_tally_mean)**2 + (uncollided_tally_err / uncollided_tally_mean)**2)

    #Hand calc at this thickness, using the unshielded flux at the detector radius (r = shield_thickness)
    mfp = mu * shield_thickness
    B_hand = float(gp_buildup(mfp))
    K_dot_uncollided_hand = flux_no_shield * np.exp(-mfp) * e_dep * MeV_per_g_s_to_Gy_per_h
    K_dot_hand = K_dot_uncollided_hand * B_hand

    print(f'\n--- Lead thickness {shield_thickness} cm ({mfp:.2f} mfp), {lead_beyond} cm of lead beyond the detector ---')
    print(f'OpenMC air kerma rate   : {K_dot_openmc * Gy_to_uGy:.4g} +/- {K_dot_openmc_err * Gy_to_uGy:.2g} uGy/h')
    print(f'Hand-calc air kerma rate: {K_dot_hand * Gy_to_uGy:.4g} uGy/h')
    print(f'Difference               : {K_dot_openmc / K_dot_hand - 1:.1%}')
    print(f'OpenMC buildup factor    : {B_openmc:.2f} +/- {B_openmc_err:.2f}')
    print(f'Hand-calc buildup factor : {B_hand:.2f}')

    return {
        'thickness_cm': shield_thickness,
        'mfp': mfp,
        'lead_beyond_cm': lead_beyond,
        'K_total_openmc_uGy_h': K_dot_openmc * Gy_to_uGy,
        'K_total_openmc_err_uGy_h': K_dot_openmc_err * Gy_to_uGy,
        'K_total_hand_uGy_h': K_dot_hand * Gy_to_uGy,
        'B_openmc': B_openmc,
        'B_openmc_err': B_openmc_err,
        'B_hand': B_hand,
    }


results_dir = Path(__file__).resolve().parent / 'Results'
results_dir.mkdir(exist_ok=True)

#Depth sweep
print('\n========== Depth sweep ==========')
sweep = pd.DataFrame([run_infinite(t, lead_beyond_detector) for t in shield_thickness])
sweep.to_csv(results_dir / 'openmc_infinite_sweep.csv', index=False)
print(f'\nSweep results saved to {results_dir / "openmc_infinite_sweep.csv"}')