import numpy as np
from matplotlib import pyplot as plt
import os, sys

############# PARAMETERS YOU CAN MODIFY TO SUITE YOUR CASES#########################

results_file_name = 'bubbles_and_stresses_selected.txt'
file_bubbles = os.path.join(os.getcwd(), results_file_name)

# Set hydrostatic pressure
p_hydro = 0.16  # Set to your actual hydrostatic pressure if known

# Set the bubble pressure
pint = 1

############# END OF THE PARAMETERS YOU CAN MODIFY       #########################



# Read the bubble data file
with open(file_bubbles, 'r') as f:
    lines = f.readlines()
    # Fix header parsing - remove empty strings from trailing comma
    header = [h.strip() for h in lines[0].strip().split(',') if h.strip()]
    data = []
    for line in lines[1:]:
        # Fix data parsing - remove empty strings
        values = [v.strip() for v in line.strip().split(',') if v.strip()]
        if values:
            data.append([float(v) for v in values])

# Check if 'Stress' column exists
if 'Stress' not in header:
    print(f"Error: 'Stress' column not found in file {file_bubbles}")
    exit()

print(f"Loaded {len(data)} bubbles")
print(f"Columns: {header}")

# Add p_hydro to 'Stress' column
for row in data:
    row[header.index('Stress')] += p_hydro

# Sort data by 'Stress' column (descending)
data.sort(key=lambda x: x[header.index('Stress')], reverse=True)

# Calculate reciprocal_stress and cdf
reciprocal_stress = np.array([(pint-p_hydro) / row[header.index('Stress')] for row in data])
cdf = np.arange(1, len(reciprocal_stress)+1) / len(reciprocal_stress)

# Save sorted data to a text file
with open('df_sorted_output.csv', 'w') as f:
    f.write(','.join(header + ['reciprocal_stress', 'cdf']) + '\n')
    for i, row in enumerate(data):
        f.write(','.join(map(str, row + [reciprocal_stress[i], cdf[i]])) + '\n')

print("Sorted data saved to 'df_sorted_output.csv'")

# Calculate bubble pressure and FGR for plotting
sig_rupt = 1.5e8  # Rupture stress in Pa
bubble_pressure = reciprocal_stress * sig_rupt
fgr_percent = 100 * cdf

# Find stress at 10% FGR
target_fgr = 10.0  # Target FGR percentage
idx_10 = np.argmin(np.abs(fgr_percent - target_fgr))
stress_at_10 = data[idx_10][header.index('Stress')]
pressure_at_10 = bubble_pressure[idx_10]
fgr_at_10 = fgr_percent[idx_10]

print(f"\n=== Stress at ~10% FGR ===")
print(f"FGR: {fgr_at_10:.2f}%")
print(f"Stress: {stress_at_10:.6f}")
print(f"Bubble pressure: {pressure_at_10:.2e} Pa")
print(f"Reciprocal stress: {reciprocal_stress[idx_10]:.6f}")

# Create plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Normalized pressure vs CDF
ax1.plot(reciprocal_stress, cdf, '-o', markersize=4, label='Bubble data')
ax1.set_xlabel('Pressure/$\\sigma_I$', fontsize=14)
ax1.set_ylabel('Cumulative bubble fraction', fontsize=14)
ax1.set_xlim(0, 3)
ax1.set_ylim(0, 1)
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_title('Normalized Stress Distribution')

# Plot 2: Bubble pressure vs FGR
ax2.plot(bubble_pressure, fgr_percent, '-o', markersize=4, label='Simulation')
# Mark the 10% FGR point
ax2.plot(pressure_at_10, fgr_at_10, 'b*', markersize=15, label=f'~10% FGR', zorder=5)
ax2.axhline(y=10, color='b', linestyle='--', alpha=0.5, linewidth=1)
ax2.axvline(x=pressure_at_10, color='b', linestyle='--', alpha=0.5, linewidth=1)
ax2.set_xlabel('Bubble pressure, Pa', fontsize=14)
ax2.set_ylabel('FGR, %', fontsize=14)
ax2.set_xlim(0, 4e8)
ax2.set_ylim(0, 100)
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.set_title('Fission Gas Release')

# Add experimental data
exp_pressure = np.array([1.04E+08, 1.08E+08, 1.11E+08, 1.14E+08, 1.17E+08, 
                          1.19E+08, 1.23E+08, 1.25E+08, 1.28E+08, 1.30E+08, 1.33E+08])
exp_fgr = np.array([0, 0.313, 0.705, 1.331, 2.427, 3.758, 5.638, 7.595, 11.823, 16.756, 22.864])
ax2.plot(exp_pressure, exp_fgr, 'r*-', markersize=10, label='Experimental')

plt.tight_layout()
plt.savefig('bubble_analysis.png')
print("Plot saved to 'bubble_analysis.png'")
plt.show()
