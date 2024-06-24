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

# Prepare data to write to file
data = np.column_stack((heights, u, v))

# Save data to file
filename = "input_ReTau2000_loglaw.txt"
np.savetxt(filename, data, header="y u v", fmt="%.8e", delimiter='\t')

print(f"Data saved to {filename}.")

