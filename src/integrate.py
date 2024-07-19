import numpy as np
from scipy.integrate import quad

# Given parameters
u_tau = 0.395  # friction velocity
nu = 0.001  # kinematic viscosity
H = 1.0  # channel height
kappa = 0.41  # von Kármán constant

# Roughness length for smooth wall
y0 = nu / u_tau

# Velocity profile function
def velocity_profile(y):
    return (u_tau / kappa) * np.log(y / y0)

# Integrate velocity profile over channel height
U_b, _ = quad(velocity_profile, y0, H)

# Normalize by channel height
U_b /= H

print(f"Bulk velocity U_b = {U_b:.4f} m/s")

