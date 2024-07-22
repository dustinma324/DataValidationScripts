import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from scipy.integrate import simpson
from scipy.integrate import simps


# Define the file path
#file_path_dns = "../chan2000/LM_Channel_2000_mean_prof.dat.txt"
#header_lines_to_skip = 71
#file_path_dns = "../chan180/MKM_Channel_180_mean_prof.txt"
file_path_dns = "../chan395/MKM_Channel_395_mean_prof.txt"
header_lines_to_skip = 24

file_path_mean = "./mean_ReTau395.dat"
#file_path_mean = "./mean_dpdpx0.08_v2.dat"

#DNS_u_tau = 4.58794e-02
#Re_tau = 1994.756
Re_tau = 395

# Smooth log-law profile parameters
smooth_Re_tau = Re_tau
smooth_mu = 0.001
smooth_delta = 1.0
smooth_u_tau = smooth_Re_tau * smooth_mu / smooth_delta
print(f"smooth_u_tau = {smooth_u_tau}")

# Compute the smooth log-law profile
def smooth_log_law_profile(y_plus, u_tau, nu):
    kappa = 0.41
    B = 5.2
    return (1/kappa) * np.log(y_plus) + B

smooth_y_plus = np.linspace(1e-1, smooth_delta, 500) * smooth_u_tau / smooth_mu
smooth_u_plus = smooth_log_law_profile(smooth_y_plus, smooth_u_tau, smooth_mu)

def read_dns_data(filename):
   columnsDNS = ["y/delta", "y^+", "U", "dU/dy", "W", "P"]
   data_DNS = pd.read_csv(file_path_dns, delimiter='\s+', skiprows=header_lines_to_skip, header=None, names=columnsDNS)
   y_plus_dns = data_DNS.iloc[1:, data_DNS.columns.get_loc("y^+")].reset_index(drop=True).astype(float)
   u_dns_plus = data_DNS.iloc[1:, data_DNS.columns.get_loc("U")].reset_index(drop=True).astype(float)
   y_over_delta = data_DNS.iloc[1:, data_DNS.columns.get_loc("y/delta")].reset_index(drop=True).astype(float)
   return y_plus_dns, u_dns_plus, y_over_delta

dns_y_plus, dns_u_plus, y_over_delta = read_dns_data(file_path_dns)
u_dns_data = dns_u_plus

# Data output for Sponging and Sounding
v_mean_data = np.zeros_like(u_dns_data)
mixing_data = np.zeros_like(u_dns_data)
potential_temp = np.full_like(u_dns_data, 300.0)
output_df_sponge = pd.DataFrame({
    'y_over_delta': y_over_delta,
    'u_mean': u_dns_data,
    'v_mean': v_mean_data
})
output_df_sounding = pd.DataFrame({
    'y_over_delta': y_over_delta,
    'potential_temp': potential_temp,
    'mixing_data':mixing_data,
    'u_mean': u_dns_data,
    'v_mean': v_mean_data
})

def read_mean_data(filename):
    columns = ["time", "height", "u", "v", "w", "rho", "theta", "tke", "col9", "col10", "col11", "col12", "col13", "col14", "col15"]
    data_toc = pd.read_csv(filename, sep='\s+', header=None, names=columns)
    toc_z = data_toc.iloc[0:, data_toc.columns.get_loc("height")].reset_index(drop=True).astype(float)
    toc_u_mean = data_toc.iloc[0:, data_toc.columns.get_loc("u")].reset_index(drop=True).astype(float)
    toc_rho = data_toc.iloc[0:, data_toc.columns.get_loc("rho")].reset_index(drop=True).astype(float)
    return toc_z, toc_u_mean, toc_rho

toc_z, toc_u_mean, toc_rho = read_mean_data(file_path_mean)

norm_by_smooth = 1

if norm_by_smooth == 1:
  # Normalizing by analytical data
  toc_y_plus = toc_z*smooth_u_tau/smooth_mu
  toc_u_plus = toc_u_mean/smooth_u_tau 
elif norm_by_smooth == 0:
  # Normalizing by simulation data
  toc_utau = np.sqrt( (smooth_mu/toc_rho[0]) * (toc_u_mean[0] / toc_z[0]) )
  toc_y_plus = toc_z * toc_utau / (smooth_mu/toc_rho[0])
  toc_u_plus = toc_u_mean / toc_utau
  print(f"toc_utau = {toc_utau}")

def compute_mean_velocity(heights, velocities):
    heights = np.array(heights)
    velocities = np.array(velocities)
    #integrated_velocity = simpson(y=velocities, x=heights)
    integrated_velocity = simps(y=velocities, x=heights)
    total_height_range = heights[-1] - heights[0]
    mean_velocity = integrated_velocity / total_height_range
    return mean_velocity

toc_ubulk = compute_mean_velocity(toc_z, toc_u_mean)
print(f"toc_bulk_u = {toc_ubulk}")

def print_output(z, toc_u_mean):
    print("Height (z):")
    print(z)
    print("\nU mean (toc_u_mean):")
    print(toc_u_mean)

#print_output(toc_z, toc_u_mean)

# Read the file into a pandas DataFrame, skipping the header lines and setting column names
try:
    output_df_sponge.to_csv('input_ReTau2000DNS_sponge.txt', sep=' ', index=False)
    output_df_sounding.to_csv('input_ReTau2000DNS_sounding.txt', sep=' ', index=False)

    # Plotting figure 1: Semilogx plot
    plt.figure(figsize=(12, 8))
    # Smooth log-law profile
    plt.semilogx(smooth_y_plus, smooth_u_plus, linestyle='--', color='r', label='Smooth Analytical Profile')
    plt.semilogx(dns_y_plus, dns_u_plus, marker='o', linestyle='-', color='b', label="LM2015 $Re_{\\tau}=$"+f"{Re_tau}")
    plt.semilogx(toc_y_plus, toc_u_plus , marker='s', linestyle='-', color='k', label="ERF Periodic $Re_{\\tau}=$"+f"{Re_tau}")
    
    # Labels
    plt.xlabel(r'$y^+$', fontsize=14, fontname='serif')
    plt.ylabel(r'$U^+$', fontsize=14, fontname='serif')
    
    # Turn off grid
    plt.grid(False)
    
    # Legend
    plt.legend(fontsize=12)
    
    # Set font to serif for all text elements
    plt.rcParams['font.family'] = 'serif'
    
    # Set global font size
    plt.rcParams.update({'font.size': 14})
    
    # Tight layout for Figure 1
    plt.tight_layout()
    
    # Save figure 1 with 600 dpi
    plt.savefig('LogLaw.png', dpi=600)
    
    # Plotting figure 2: Second plot
    plt.figure(figsize=(12, 8))
    plt.plot(u_dns_data, y_over_delta, marker='o', linestyle='-', color='b', label="LM2015 $Re_{\\tau}=$"+f"{Re_tau}")
    plt.plot(toc_u_mean/toc_ubulk, toc_z, marker='s', linestyle='-', color='k', label="ERF Periodic $Re_{\\tau}=$"+f"{Re_tau}")
    
    # Labels
    plt.xlabel(r'$\overline{U}/U_{B}$', fontsize=14, fontname='serif')
    plt.ylabel(r'$y/\delta$', fontsize=14, fontname='serif')
    
    # Set x and y limits
    plt.xlim(0.0, 1.4)
    plt.ylim(0.0, 1.0)
    
    # Legend
    plt.legend(fontsize=12)
    
    # Set font to serif for all text elements
    plt.rcParams['font.family'] = 'serif'
    
    # Set global font size
    plt.rcParams.update({'font.size': 14})
    
    # Tight layout for Figure 2
    plt.tight_layout()
    
    # Save figure 2 with 600 dpi
    plt.savefig('MeanVelocity.png', dpi=600)
    
    # Show plots
    plt.show()
    
except FileNotFoundError:
    print(f"The file at {file_path} does not exist.")
except pd.errors.EmptyDataError:
    print("The file is empty.")
except pd.errors.ParserError:
    print("Error parsing the file.")

