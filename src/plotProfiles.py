import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define the file path
file_path_dns = "../chan2000/LM_Channel_2000_mean_prof.dat.txt"
file_path_mean = "./mean.dat"

# Number of header rows to skip
header_lines_to_skip = 71

# Column names
column_names = ["y/delta", "y^+", "U", "dU/dy", "W", "P"]

# Additional variables
nx = 4096
ny = 768
nz = 3072
Lx = 8 * np.pi
Lz = 3 * np.pi
n = 7
nu = 2.30000e-05
delta = 1.000
U_mean = 1.000
u_tau = 4.58794e-02
Re_tau = 1994.756

# Smooth log-law profile parameters
smooth_Re_tau = 2000
smooth_u_tau = 0.2
smooth_nu = 0.0001
smooth_delta = 1.0

# Compute the smooth log-law profile
def smooth_log_law_profile(y_plus, u_tau, nu):
    kappa = 0.41
    B = 5.2
    return (1/kappa) * np.log(y_plus) + B

smooth_y_plus = np.linspace(1e-3, smooth_delta, 500) * smooth_u_tau / smooth_nu
smooth_u_plus = smooth_log_law_profile(smooth_y_plus, smooth_u_tau, smooth_nu)

def read_dns_data(filename):
   columnsDNS = ["y/delta", "y^+", "U", "dU/dy", "W", "P"]
   data_DNS = pd.read_csv(file_path_dns, delimiter='\s+', skiprows=header_lines_to_skip, header=None, names=column_names)
   y_plus = data_DNS.iloc[1:, data_DNS.columns.get_loc("y^+")].reset_index(drop=True).astype(float)
   u_plus = data_DNS.iloc[1:, data_DNS.columns.get_loc("U")].reset_index(drop=True).astype(float)
   y_over_delta = data_DNS.iloc[1:, data_DNS.columns.get_loc("y/delta")].reset_index(drop=True).astype(float)
   return y_plus, u_plus, y_over_delta

dns_y_plus, dns_u_plus, y_over_delta = read_dns_data(file_path_dns)
u_dns_data = dns_u_plus * u_tau
v_mean_data = np.zeros_like(u_dns_data)
output_df = pd.DataFrame({
    'y_over_delta': y_over_delta,
    'u_mean': u_dns_data,
    'v_mean': v_mean_data
})

def read_mean_data(filename):
    columns = ["time", "height", "u", "v", "w", "rho", "theta", "tke", "col9", "col10", "col11", "col12", "col13", "col14"]
    data_toc = pd.read_csv(filename, sep='\s+', header=None, names=columns)
    toc_z = data_toc.iloc[1:, data_toc.columns.get_loc("height")].reset_index(drop=True).astype(float)
    toc_u_mean = data_toc.iloc[1:, data_toc.columns.get_loc("u")].reset_index(drop=True).astype(float)
    return toc_z, toc_u_mean

toc_z, toc_u_mean = read_mean_data(file_path_mean)
toc_y_plus = toc_z*smooth_u_tau/smooth_nu
toc_u_plus = toc_u_mean/smooth_u_tau

def print_output(z, toc_u_mean):
    print("Height (z):")
    print(z)
    print("\nU mean (toc_u_mean):")
    print(toc_u_mean)

#print_output(toc_z, toc_u_mean)

# Read the file into a pandas DataFrame, skipping the header lines and setting column names
try:
    output_df.to_csv('input_ReTau2000DNS.txt', sep=' ', index=False)

    # Plotting figure 1: Semilogx plot
    plt.figure(figsize=(12, 8))
    # Smooth log-law profile
    plt.semilogx(smooth_y_plus, smooth_u_plus, linestyle='--', color='r', label='Smooth Analytical Profile')

    # DNS profile
    plt.semilogx(dns_y_plus, dns_u_plus, marker='o', linestyle='-', color='b', label="LM2015 $Re_{\\tau}=2000$")

    # Turbulent Open Channel Periodic Data
    plt.semilogx(toc_y_plus, toc_u_plus , marker='s', linestyle='-', color='k', label="ERF Periodic")
    
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
    plt.plot(u_dns_data, y_over_delta, marker='o', linestyle='-', color='b', label="LM2015 $Re_{\\tau}=2000$")
    plt.plot(toc_u_mean, toc_z, marker='s', linestyle='-', color='k', label="ERF Periodic")
    
    # Labels
    plt.xlabel(r'$\overline{U}$', fontsize=14, fontname='serif')
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

