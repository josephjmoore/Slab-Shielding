import os
from pathlib import Path
import numpy as np
 

cross_sections = os.environ.get(
    'OPENMC_CROSS_SECTIONS', '/home/asus/nuclear_data/endfb-viii.0-hdf5/cross_sections.xml')
 
#NIST air table (energy MeV, mu/rho cm^2/g, mu_en/rho cm^2/g)
air_data_file = Path(__file__).resolve().parent.parent.parent / 'Data' / 'TargetAirMu.txt'
 

#Unit conversion factors

MeV_to_J = 1.602176634e-13 #J per MeV
g_to_kg = 1000.0 #g per kg (converts per g to per kg)
s_to_h = 3600.0 #s per h (converts per s to per h)
Gy_to_uGy = 1e6 #uGy per Gy
MeV_to_eV = 1e6 #eV per MeV (OpenMC uses eV)
 
#Source
A = 1e9 #1GBq
E = 0.662 #MeV
#Number of photons emitted per decay
n_photons = 0.851
#Photon emission rate in photons/s
S = A * n_photons


#Geometry (cm)
r_detector = 100.0 #cm
r_in = 10.0 #cm
r_det_in = r_detector - 0.5 #cm
r_det_out = r_detector + 0.5 #cm
r_end = r_det_out + 10.0 #cm
#Volume of the detector shell in cm^3 to turn the flux tally into flux
V_det = (4.0 / 3.0) * np.pi * (r_det_out**3 - r_det_in**3)
 
#Flux with no shield at the detector in photons/cm^2/s
flux_no_shield = S / (4 * np.pi * r_detector**2)
 

#Lead shield
lead_density = 11.34 #g/cm^3
#Mass attenuation coefficient from NIST XCOM at 0.662 MeV
#(Total Attenuation WITHOUT Coherent Scattering, to match the buildup tables)
mu_rho = 0.1035 #cm^2/g
#Linear attenuation coefficient in cm^-1
mu = mu_rho * lead_density #cm^-1
 

#Air 'detector' response function
def loglog_interp(x, xp, fp):
    """Log-log interpolation of fp(xp) at x (xp must be increasing)."""
    return np.exp(np.interp(np.log(x), np.log(xp), np.log(fp)))
 
#Whole NIST air table: energies in MeV and mu_en/rho in cm^2/g
air_energy, air_mu_en_rho = np.loadtxt(air_data_file, usecols=(0, 2), unpack=True)
 
#mu_en/rho for air at the source energy in cm^2/g
mu_en_rho = loglog_interp(E, air_energy, air_mu_en_rho)
 
#Energy deposited per unit flux in air 'detector' in MeV*cm^2/g
e_dep = mu_en_rho * E
 
#Air kerma rate with no shield in Gy/h (should be ~75.7 uGy/h)
K_dot_no_shield = flux_no_shield * e_dep * MeV_to_J * g_to_kg * s_to_h
 
#GP buildup factor for lead, exposure [Shultis and Faw, 2000]
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
#Fraction of the way from 0.6 to 0.8 MeV
f = (E - 0.6) / (0.8 - 0.6)
#Parameters interpolated to the source energy
b = b_06 + f * (b_08 - b_06)
c = c_06 + f * (c_08 - c_06)
a = a_06 + f * (a_08 - a_06)
Xk = Xk_06 + f * (Xk_08 - Xk_06)
d = d_06 + f * (d_08 - d_06)
 
 
def gp_buildup(x, b=b, c=c, a=a, Xk=Xk, d=d):
    """GP buildup factor at depth x (mean free paths).
 
    x can be a single value or a numpy array of depths.
    Uses the parameters interpolated to the source energy unless others are given.
    """
    x = np.asarray(x, dtype=float)
    tanh_m2 = np.tanh(-2.0)
    K = c * x**a + d * (np.tanh(x / Xk - 2.0) - tanh_m2) / (1.0 - tanh_m2)
    #Where K is 1 the GP formula divides by zero, so use the linear form there
    K_is_one = np.isclose(K, 1.0)
    K_safe = np.where(K_is_one, 2.0, K)
    B_gp = 1.0 + (b - 1.0) * (K_safe**x - 1.0) / (K_safe - 1.0)
    B_linear = 1.0 + (b - 1.0) * x
    return np.where(K_is_one, B_linear, B_gp)
