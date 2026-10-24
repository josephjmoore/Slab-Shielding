from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from params import mu, K_dot_no_shield, gp_buildup, Gy_to_uGy

csv_path = Path(__file__).resolve().parent / "Results" / "openmc_sweep.csv"
csv = pd.read_csv(csv_path)

csv_inf_path = Path(__file__).resolve().parent / "Results" / "openmc_infinite_sweep.csv"
csv_inf = pd.read_csv(csv_inf_path)

shield_thickness = csv["thickness_cm"]
K_total_openmc_uGy_h = csv["K_total_openmc_uGy_h"]
K_total_openmc_err_uGy_h = csv["K_total_openmc_err_uGy_h"]
K_total_hand_uGy_h = csv["K_total_hand_uGy_h"]
K_uncollided_openmc_uGy_h = csv["K_uncollided_openmc_uGy_h"]
K_uncollided_openmc_err_uGy_h = csv["K_uncollided_openmc_err_uGy_h"]
K_uncollided_hand_uGy_h = csv["K_uncollided_hand_uGy_h"]
mfp = csv["mfp"]
B_openmc = csv["B_openmc"]
B_openmc_err = csv["B_openmc_err"]

#Hand calculation as a smooth curve (it's a formula, so evaluate it on a fine grid)
t_fine = np.linspace(0, shield_thickness.max(), 200) #cm
mfp_fine = mu * t_fine
B_hand_fine = gp_buildup(mfp_fine)
K_hand_fine = K_dot_no_shield * np.exp(-mfp_fine) * B_hand_fine * Gy_to_uGy

#Ratio of OpenMC to hand calculation at each simulated thickness
ratio_total = K_total_openmc_uGy_h / K_total_hand_uGy_h
ratio_total_err = K_total_openmc_err_uGy_h / K_total_hand_uGy_h
ratio_uncollided = K_uncollided_openmc_uGy_h / K_uncollided_hand_uGy_h
ratio_uncollided_err = K_uncollided_openmc_err_uGy_h / K_uncollided_hand_uGy_h

#Effective infinite thickness for build up comparison - validates openmc
mfp_inf = csv_inf["mfp"]
B_openmc_inf = csv_inf["B_openmc"]
B_openmc_err_inf = csv_inf["B_openmc_err"]

#Colours
colour_openmc = "#83d62a"      
colour_hand = '#eb6834' 
colour_total = "#ff00fb" 
colour_uncollided = '#3d3d3a' 
colour_reference = '#8a8a85' 
colour_grid = '#d9d9d6'
colour_text = '#3d3d3a'

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.5))

#Left: air kerma rate against shield thickness
ax1.plot(t_fine, K_hand_fine,
         color=colour_hand, linewidth=2, linestyle='--', label='Hand Calculation')
ax1.errorbar(shield_thickness, K_total_openmc_uGy_h, yerr=K_total_openmc_err_uGy_h,
             fmt='x', color=colour_openmc, markersize=6, capsize=3, label='OpenMC')
ax1.set_yscale('log')
ax1.set_xlabel('Lead thickness (cm)')
ax1.set_ylabel('Air kerma rate at 1 m (µGy/h)')
ax1.set_title('Cs-137 (1 GBq) behind lead', loc='left')
ax1.legend(frameon=False)

#Middle: buildup factor against depth in mean free paths
ax2.plot(mfp_fine, B_hand_fine,
         color=colour_hand, linewidth=2, linestyle='--', label='Hand Calculation (GP)')
ax2.errorbar(mfp, B_openmc, yerr=B_openmc_err,
             fmt='x', color=colour_openmc, markersize=6, capsize=3, label='OpenMC')
ax2.errorbar(mfp_inf, B_openmc_inf, yerr=B_openmc_err_inf,
             fmt='x', color=colour_text, markersize=6, capsize=3, label='OpenMC Infinite Effective Thickness')
ax2.set_xlabel('Lead thickness (mean free paths)')
ax2.set_ylabel('Exposure buildup factor B')
ax2.set_title('Buildup factor, lead, 0.662 MeV', loc='left')
ax2.set_ylim(bottom=1.0)
ax2.legend(frameon=False)

#Right: ratio of OpenMC to hand calculation
ax3.axhline(1.0, color=colour_reference, linewidth=1, linestyle='--')
ax3.errorbar(shield_thickness, ratio_uncollided, yerr=ratio_uncollided_err,
             fmt='x', color=colour_uncollided, markerfacecolor='none', markersize=6,
             capsize=3, label='Uncollided')
ax3.errorbar(shield_thickness, ratio_total, yerr=ratio_total_err,
             fmt='x', color=colour_total, markersize=6, capsize=3, label='Total')
ax3.set_xlabel('Lead thickness (cm)')
ax3.set_ylabel('OpenMC / Hand Calculation')
ax3.set_title('OpenMC relative to hand calculation', loc='left')
ax3.legend(frameon=False)

#Keep gridlines and axes recessive so the data stands out
for ax in (ax1, ax2, ax3):
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

#Save the figure in the Results folder next to this script
fig.savefig(csv_path.parent / 'openmc_vs_hand_calc.png', dpi=200)
print(f'Plot saved to {csv_path.parent / "openmc_vs_hand_calc.png"}')

plt.show()