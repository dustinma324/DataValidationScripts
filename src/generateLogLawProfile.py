import numpy as np

# Given parameters
Re_tau = 2000
height_start = 1e-2
height_end = 10.0
num_points = 33
kappa = 0.41  # Von Kármán constant
B = 5.2

# Friction velocity (u_tau)
u_tau = 1.0  # Set to 1 for simplicity, actual value would depend on Re_tau

# Generate heights (y)
heights = np.linspace(height_start, height_end, num_points)

# Calculate u values using smooth wall log law of the wall formula
u = u_tau * (1./ kappa * np.log(heights * Re_tau) + B) 

# Create v values (set all v values to 0)
v = np.zeros_like(heights)
mixing_ratio = np.zeros_like(heights)
potential_temp = np.full_like(heights, 300.0)

# Prepare data to write to file
data_sounding = np.column_stack((heights, potential_temp, mixing_ratio, u, v))
data_sponge = np.column_stack((heights, u, v))

# Save data to file
filename = "input_ReTau2000Ana"
filetype = ".txt"

np.savetxt(filename+"_sounding"+filetype, data_sounding, header="y theta mr u v", fmt="%.8e", delimiter='\t')
print(f"Data saved to {filename}_sounding{filetype}.")
np.savetxt(filename+"_sponge"+filetype, data_sponge, header="y u v", fmt="%.8e", delimiter='\t')
print(f"Data saved to {filename}_sponge{filetype}.")

