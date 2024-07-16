import numpy as np
import matplotlib.pyplot as plt

# Given parameters
Re_tau = 2000
nu = 0.0001  # Kinematic viscosity in m^2/s
delta = 8.0  # Half channel height in meters

# Calculating friction velocity
u_tau = (Re_tau * nu) / delta

# Constants for the log-law
kappa = 0.41
B = 5.2

# Heights from 0 to 8 meters
y = np.linspace(0.001, 8.0, 500)  # Avoid y=0 to prevent log(0) error
y_plus = (y * u_tau) / nu

# Calculate U^+ using the log-law formula
U_plus = (1 / kappa) * np.log(y_plus) + B

# Convert U^+ to U
U = U_plus * u_tau

# Plotting the velocity profile
plt.figure(figsize=(10, 6))
plt.plot(U, y, label=r'Smooth Log-Law Profile')
plt.xlabel(r'$U$ (m/s)', fontsize=14, fontname='serif')
plt.ylabel(r'$y$ (m)', fontsize=14, fontname='serif')
plt.title(r'Smooth Log-Law Velocity Profile for $Re_{\tau}=2000$', fontsize=16, fontname='serif')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

