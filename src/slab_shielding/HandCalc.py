from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from params import mu, K_dot_no_shield, Gy_to_uGy, gp_buildup, b, c, a, Xk, d

#Shield thicknesses to sweep in cm
shield_thickness = np.linspace(0.0, 10.0, 41) #cm

#Number of mean free paths for each thickness
mfp_s = mu * shield_thickness

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

#Plots
#Colours
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


