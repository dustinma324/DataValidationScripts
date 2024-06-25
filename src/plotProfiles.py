import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define the file path
file_path = "../chan2000/LM_Channel_2000_mean_prof.dat.txt"

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

y_plus_smooth = np.linspace(1e-3, smooth_delta, 500) * smooth_u_tau / smooth_nu
u_plus_smooth = smooth_log_law_profile(y_plus_smooth, smooth_u_tau, smooth_nu)

# Read the file into a pandas DataFrame, skipping the header lines and setting column names
try:
    df = pd.read_csv(file_path, delimiter='\s+', skiprows=header_lines_to_skip, header=None, names=column_names)
    
    # Extracting the required columns and converting to floats
    y_plus = df.iloc[1:, df.columns.get_loc("y^+")].reset_index(drop=True).astype(float)
    u_plus = df.iloc[1:, df.columns.get_loc("U")].reset_index(drop=True).astype(float)
    y_over_delta = df.iloc[1:, df.columns.get_loc("y/delta")].reset_index(drop=True).astype(float)
    
    # Calculate u_mean (u_mean_data)
    u_mean_data = u_plus * u_tau
    
    # Calculate v_mean (all zeros)
    v_mean_data = np.zeros_like(u_mean_data)
    
    # Create a DataFrame with desired columns
    output_df = pd.DataFrame({
        'y_over_delta': y_over_delta,
        'u_mean': u_mean_data,
        'v_mean': v_mean_data
    })
    
    # Output DataFrame to TXT file with tab-separated values
    output_df.to_csv('input_ReTau2000.txt', sep='\t', index=False)
    
    # Plotting figure 1: Semilogx plot
    plt.figure(figsize=(12, 8))
    # Smooth log-law profile
    plt.semilogx(y_plus_smooth, u_plus_smooth, linestyle='--', color='r', label='Smooth Analytical Profile')

    # DNS profile
    plt.semilogx(y_plus, u_plus, marker='o', linestyle='-', color='b', label="LM2015 $Re_{\\tau}=2000$")
    
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
    plt.plot(u_mean_data, y_over_delta, marker='o', linestyle='-', color='g', label="LM2015 $Re_{\\tau}=2000$")
    
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

