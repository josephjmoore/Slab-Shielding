import openmc
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

#Create Shield Material
shield = openmc.Material(name='lead shield')
shield.add_element('Pb', 1.0)
shield.set_density('g/cm3', 11.34)
shield_thickness = 5 #cm

#Source to detector distance in cm
r_detector = 100.0 #cm
#Inner radius of the lead shell in cm
r_in = 10.0 #cm


def build_geometry(shield_thickness):
    """Build the spherical geometry for one lead thickness (cm).

    Returns the geometry and the detector cell (needed for the tallies).
    """

    r_out = r_in + shield_thickness 
    r_det_in = r_detector - 0.5 
    r_det_out = r_detector + 0.5
    r_end = r_det_out + 10.0

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
geometry.export_to_xml()