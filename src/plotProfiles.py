import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

# Read the file into a pandas DataFrame, skipping the header lines and setting column names
try:
    df = pd.read_csv(file_path, delimiter='\s+', skiprows=header_lines_to_skip, header=None, names=column_names)
    
    # Extracting the required columns and converting to floats
    y_plus = df.iloc[1:, df.columns.get_loc("y^+")].reset_index(drop=True).astype(float)
    u_plus = df.iloc[1:, df.columns.get_loc("U")].reset_index(drop=True).astype(float)
    y_over_delta = df.iloc[1:, df.columns.get_loc("y/delta")].reset_index(drop=True).astype(float)
    
    # Calculate u_mean
    u_mean = u_plus * u_tau
    
    # Plotting figure 1: Semilogx plot
    plt.figure(figsize=(12, 8))
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
    plt.savefig('LogLaw.png', dpi=600)
    
    # Plotting figure 2: Second plot
    plt.figure(figsize=(12, 8))
    plt.plot(u_mean, y_over_delta, marker='o', linestyle='-', color='g', label="LM2015 $Re_{\\tau}=2000$")
    
    # Labels
    plt.xlabel(r'$u_{mean}$', fontsize=14, fontname='serif')
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
    plt.savefig('MeanVelocity.png', dpi=600)
    
    # Show plots
    plt.show()
    
except FileNotFoundError:
    print(f"The file at {file_path} does not exist.")
except pd.errors.EmptyDataError:
    print("The file is empty.")
except pd.errors.ParserError:
    print("Error parsing the file.")

