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
r_detector = 100.0 #cm

#Flux no shield in photons/cm^2/s
flux_no_shield = S / (4 * np.pi * r_detector**2)

#Create Shield Material
shield = openmc.Material(name='lead shield')
shield.add_element('Pb', 1.0)
shield.set_density('g/cm3', 11.34)
#Shield thicknesses to sweep in cm
shield_thickness = np.linspace(0.0, 10.0, 41) #cm
#Attenuation coefficient from NIST XCOM database for 0.662 MeV photons in lead
mu_rho = 0.1035 #cm^2/g
mu = mu_rho * shield.density #cm^-1

#Number of mean free paths for each thickness
mfp_s = mu * shield_thickness

#Target Air Attenuation Coefficients for 'Detector' Response Function
target_air_mu = np.loadtxt(
    Path(__file__).resolve().parent.parent.parent / 'Data' / 'TargetAirMu.txt',
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
    """GP buildup factor at depth x (mean free paths).

    x can be a single value or a numpy array of depths.
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

#Energy deposited per unit flux in air 'detector' in MeV*cm^2/g
e_dep = mu_en_rho * E

#Air kerma rate with no shield in Gy/h
K_dot_no_shield = flux_no_shield * e_dep * MeV_to_J * g_to_kg * s_to_h

#Buildup factor for each thickness
B = gp_buildup(mfp_s, b, c, a, Xk, d)

#Air kerma rate from uncollided photons only, in Gy/h
K_dot_uncollided = K_dot_no_shield * np.exp(-mfp_s)
#Total air kerma rate (uncollided + scattered), in Gy/h
K_dot_total = K_dot_uncollided * B

#Print every 1 cm (every 4th point)
print('\nThickness (cm)   mfp     B      Uncollided (uGy/h)   Total (uGy/h)')
for i in range(0, len(shield_thickness), 4):
    print(f'{shield_thickness[i]:10.1f} {mfp_s[i]:9.2f} {B[i]:7.3f}'
          f' {K_dot_uncollided[i] * Gy_to_uGy:16.4e} {K_dot_total[i] * Gy_to_uGy:16.4e}')

#----------------------------------------------------------------------------
#Plots
#----------------------------------------------------------------------------
#Colours (colour-blind safe pair)
colour_total = '#2a78d6'      #blue
colour_uncollided = '#eb6834' #orange
colour_grid = '#d9d9d6'
colour_text = '#3d3d3a'

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

#Left: air kerma rate against shield thickness
ax1.plot(shield_thickness, K_dot_total * Gy_to_uGy,
         color=colour_total, linewidth=2, label='Total (with buildup)')
ax1.plot(shield_thickness, K_dot_uncollided * Gy_to_uGy,
         color=colour_uncollided, linewidth=2, linestyle='--', label='Uncollided only')
ax1.set_yscale('log')
ax1.set_xlabel('Lead thickness (cm)')
ax1.set_ylabel('Air kerma rate at 1 m (µGy/h)')
ax1.set_title('Cs-137 (1 GBq) behind lead', loc='left')
ax1.legend(frameon=False)

#Right: buildup factor against depth in mean free paths
ax2.plot(mfp_s, B, color=colour_total, linewidth=2)
ax2.set_xlabel('Lead thickness (mean free paths)')
ax2.set_ylabel('Exposure buildup factor B')
ax2.set_title('GP buildup factor, lead, 0.662 MeV', loc='left')
ax2.set_ylim(bottom=1.0)

#Keep gridlines and axes recessive so the data stands out
for ax in (ax1, ax2):
    ax.grid(True, which='major', color=colour_grid, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(colour_grid)
    ax.tick_params(colors=colour_text)
    ax.xaxis.label.set_color(colour_text)
    ax.yaxis.label.set_color(colour_text)
    ax.title.set_color(colour_text)

fig.tight_layout()

#Save the figure in a Results folder next to this script
results_dir = Path(__file__).resolve().parent / 'Results'
results_dir.mkdir(exist_ok=True)
fig.savefig(results_dir / 'hand_calc_sweep.png', dpi=200)
print(f'\nPlot saved to {results_dir / "hand_calc_sweep.png"}')

plt.show()




