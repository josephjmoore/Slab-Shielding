import openmc
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from params import (cross_sections, S, E, MeV_to_eV, flux_no_shield, lead_density,
                    r_detector, r_in, r_det_in, r_det_out, r_end, V_det, mu, air_energy, 
                    air_mu_en_rho, MeV_to_J, g_to_kg, s_to_h, e_dep, K_dot_no_shield, gp_buildup, Gy_to_uGy)

#Create Shield Material
shield = openmc.Material(name='lead shield')
shield.add_element('Pb', 1.0)
shield.set_density('g/cm3', lead_density)
shield_thickness = 5.0 #cm


def build_geometry(shield_thickness):
    """Build the spherical geometry for one lead thickness (cm).

    Returns the geometry and the detector cell (needed for the tallies).
    """

    r_out = r_in + shield_thickness

    #Surfaces: one sphere per radius, shared by the cells either side of it
    shield_inner_surface = openmc.Sphere(r=r_in)
    detector_inner_surface = openmc.Sphere(r=r_det_in)
    detector_outer_surface = openmc.Sphere(r=r_det_out)
    #Only the outermost surface is a vacuum boundary (kills escaping particles)
    outer_boundary = openmc.Sphere(r=r_end, boundary_type='vacuum')

    if shield_thickness > 0:
        shield_outer_surface = openmc.Sphere(r=r_out)

        source_void_cell = openmc.Cell(name='source void', region=-shield_inner_surface)

        shield_cell = openmc.Cell(name='lead shield', region=+shield_inner_surface & -shield_outer_surface)
        shield_cell.fill = shield

        gap_cell = openmc.Cell(name='gap', region=+shield_outer_surface & -detector_inner_surface)
        cells = [source_void_cell, shield_cell, gap_cell]
    else:
        #No lead
        gap_cell = openmc.Cell(name='source void and gap', region=-detector_inner_surface)
        cells = [gap_cell]

    detector_cell = openmc.Cell(name='detector', region=+detector_inner_surface & -detector_outer_surface)

    outer_void_cell = openmc.Cell(name='outer void', region=+detector_outer_surface & -outer_boundary)
    cells += [detector_cell, outer_void_cell]

    root = openmc.Universe(cells=cells)
    geometry = openmc.Geometry(root)
    return geometry, detector_cell


geometry, detector_cell = build_geometry(shield_thickness)

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
materials = openmc.Materials([shield])
openmc.config['cross_sections'] = cross_sections

model = openmc.Model(geometry=geometry, materials=materials, settings=settings, tallies=tallies)

run_dir = Path(__file__).resolve().parent / 'runs' / f't_{shield_thickness}cm'
statepoint_path = model.run(cwd=run_dir)
with openmc.StatePoint(statepoint_path) as sp:
    result = sp.get_tally(name='detector tally')
    flux_tally = result.mean.item()
    flux_tally_err = result.std_dev.item()

    uncollided_result = sp.get_tally(name='uncollided tally')
    uncollided_tally_mean = uncollided_result.mean.item()
    uncollided_tally_err = uncollided_result.std_dev.item()

    dose_result = sp.get_tally(name='air kerma tally')
    dose_tally_mean = dose_result.mean.item()
    dose_tally_err = dose_result.std_dev.item()

#Total flux (uncollided + scattered) - no hand-calc equivalent
flux_openmc = flux_tally / V_det * S
#Uncollided flux
uncollided_flux_openmc = uncollided_tally_mean / V_det * S
uncollided_flux_hand = flux_no_shield * np.exp(-mu * shield_thickness)
#Linear attenuation coefficient OpenMC effectively used
mu_openmc = -np.log(uncollided_tally_mean) / shield_thickness
#Total air kerma rate (uncollided + scattered), in Gy/h
K_dot_openmc = dose_tally_mean / V_det * S * MeV_to_J * g_to_kg * s_to_h
K_dot_openmc_err = dose_tally_err / V_det * S * MeV_to_J * g_to_kg * s_to_h
#Uncollided air kerma rate, in Gy/h (all at 662 keV, so e_dep applies)
K_dot_uncollided_openmc = uncollided_flux_openmc * e_dep * MeV_to_J * g_to_kg * s_to_h
#Monte Carlo buildup factor
B_openmc = K_dot_openmc / K_dot_uncollided_openmc

#Hand calc at this thickness
mfp = mu * shield_thickness
B_hand = gp_buildup(mfp)
K_dot_uncollided_hand = K_dot_no_shield * np.exp(-mfp)
K_dot_hand = K_dot_uncollided_hand * B_hand


print(f'OpenMC air kerma rate   : {K_dot_openmc * Gy_to_uGy:.2f} uGy/h')
print(f'OpenMC air kerma rate error : {K_dot_openmc_err * Gy_to_uGy:.2f} uGy/h')
print(f'Hand-calc air kerma rate: {K_dot_hand * Gy_to_uGy:.2f} uGy/h')
print(f'Difference               : {K_dot_openmc / K_dot_hand - 1:.1%}')
print(f'OpenMC buildup factor    : {B_openmc:.2f}')
print(f'Hand-calc buildup factor : {B_hand:.2f}')