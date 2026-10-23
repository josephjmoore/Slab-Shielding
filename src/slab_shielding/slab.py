import openmc
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

cross_sections = os.environ.get(
    'OPENMC_CROSS_SECTIONS', '/home/asus/nuclear_data/endfb-viii.0-hdf5/cross_sections.xml')

#Unit conversion factors
MeV_to_J = 1.602176634e-13 #J per MeV
g_to_kg = 1000.0 #g per kg (converts per g to per kg)
s_to_h = 3600.0 #s per h (converts per s to per h)
Gy_to_uGy = 1e6 #uGy per Gy

#Source activity in Bq
A = 1e9 #1GBq
#Photon energy in MeV
E = 0.662 #MeV
#Number of photons emitted per decay
n_photons = 0.851
#Photon emission rate in photons/s
S = A * n_photons

#Source to detector distance in cm
r = 100.0 #cm

#Flux no shield in photons/cm^2/s
flux_no_shield = S / (4 * np.pi * r**2)

#Create Shield Material
shield = openmc.Material(name='lead shield')
shield.add_element('Pb', 1.0)
#Shield density in g/cm^3
shield.set_density('g/cm3', 11.34)
#Shield thickness in cm
shield_thickness = 5.0 #cm
#Attenuation coefficient from NIST XCOM database for 0.662 MeV photons in lead
mu_rho = 0.1035 #cm^2/g
mu = mu_rho * shield.density #cm^-1

#number of mean free paths
mfp = mu*shield_thickness

#Target Air Attenuation Coefficients for 'Detector' Response Function
target_air_mu = np.loadtxt(
    Path(__file__).resolve().parent / 'Data' / 'TargetAirMu.txt',
)
# log log interpolation of attenuation coefficients for target air between 0.6 and 0.8 MeV
flog = (np.log(E/0.6))/(np.log(0.8/0.6))
mu_en_rho = target_air_mu[24, 2] * ((target_air_mu[25, 2] / target_air_mu[24, 2]))**flog

#Buildup factor parameters for lead [Shultis and Faw, 2000]
#Require interpolation between the values for 0.6 MeV and 0.8 MeV
f = (E - 0.6) / (0.8 - 0.6)
#0.6 MeV
b_06 = 1.269
c_06 = 0.685
a_06 = 0.089
Xk_06 = 13.78
d_06 = -0.0440
#0.8 MeV
b_08 = 1.329
c_08 = 0.759
a_08 = 0.065
Xk_08 = 13.69
d_08 = -0.0317
#0.662 is 31% of the way from 0.6 to 0.8 MeV
b = b_06 + f * (b_08 - b_06)
c = c_06 + f * (c_08 - c_06)
a = a_06 + f * (a_08 - a_06)
Xk = Xk_06 + f * (Xk_08 - Xk_06)
d = d_06 + f * (d_08 - d_06)

def gp_buildup(x, b, c, a, Xk, d):
    """GP buildup factor at depth x (mean free paths)."""
    tanh_m2 = np.tanh(-2.0)
    K = c * x**a + d * (np.tanh(x / Xk - 2.0) - tanh_m2) / (1.0 - tanh_m2)
    if np.isclose(K, 1.0):
        return 1.0 + (b - 1.0) * x
    return 1.0 + (b - 1.0) * (K**x - 1.0) / (K - 1.0)

#Energy deposited per unit flux in air 'detector' in MeV*cm^2/g
e_dep = mu_en_rho * E

#Air kerma rate with no shield in MeV/g/s (sanity check)
K_dot_no_shield = flux_no_shield * e_dep
#Convert MeV/g/s to Gy/h
K_dot_no_shield = K_dot_no_shield * MeV_to_J * g_to_kg * s_to_h

#Air kerma rate behind the shield in MeV/g/s
K_dot = flux_no_shield * gp_buildup(mfp, b, c, a, Xk, d) * np.exp(-mfp) * e_dep
#Convert MeV/g/s to Gy/h
K_dot = K_dot * MeV_to_J * g_to_kg * s_to_h

print(f'Air kerma rate, no shield    : {K_dot_no_shield * Gy_to_uGy:.2f} uGy/h (expect ~75.7)')
print(f'Air kerma rate, behind shield: {K_dot * Gy_to_uGy:.4f} uGy/h')