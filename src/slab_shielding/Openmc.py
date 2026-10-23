import openmc
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from params import (cross_sections, S, E, MeV_to_eV, flux_no_shield, lead_density,
                    r_detector, r_in, r_det_in, r_det_out, r_end, V_det, mu)

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

tallies = openmc.Tallies([tally, uncollided_tally])
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

#Total flux (uncollided + scattered) - no hand-calc equivalent
flux_openmc = flux_tally / V_det * S
#Uncollided flux
uncollided_flux_openmc = uncollided_tally_mean / V_det * S
uncollided_flux_hand = flux_no_shield * np.exp(-mu * shield_thickness)
#Linear attenuation coefficient OpenMC effectively used
mu_openmc = -np.log(uncollided_tally_mean) / shield_thickness

print(f'Flux tally (per source photon)      : {flux_tally:.4e} +/- {flux_tally_err:.4e} cm')
print(f'Uncollided tally (per source photon): {uncollided_tally_mean:.4e} +/- {uncollided_tally_err:.4e} cm')
print(f'Uncollided relative error           : {uncollided_tally_err / uncollided_tally_mean:.2%}')
print(f'OpenMC total flux       : {flux_openmc:.2f} photons/cm^2/s')
print(f'OpenMC uncollided flux  : {uncollided_flux_openmc:.2f} photons/cm^2/s')
print(f'Hand-calc uncollided    : {uncollided_flux_hand:.2f} photons/cm^2/s')
print(f'Difference              : {uncollided_flux_openmc / uncollided_flux_hand - 1:.1%}')
print(f'mu (OpenMC)             : {mu_openmc:.4f} cm^-1')
print(f'mu (hand calc, XCOM)    : {mu:.4f} cm^-1')