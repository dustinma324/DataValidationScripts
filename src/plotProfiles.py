import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps

# Define the file path
file_path = "/home/tma/Desktop/DataValidationScripts/chan395/profiles/chan395.means"

# Number of header rows to skip
header_lines_to_skip = 24

# Read the file into a pandas DataFrame, skipping the header lines and setting column names
try:
    df = pd.read_csv(file_path, delimiter='\s+', skiprows=header_lines_to_skip, header=None, names=["y", "y+", "Umean", "col4", "col5", "col6", "col7"])  # Use '\s+' for any whitespace delimiter
    
    # Extracting the required columns and converting to floats
    y = df["y"].iloc[1:].reset_index(drop=True).astype(float)
    Umean = df["Umean"].iloc[1:].reset_index(drop=True).astype(float)
    
    # Calculate Ubulk (normalized integral of Umean with respect to y)
    H = 1.0  # Assuming H is the height and normalizing factor
    Ubulk = simps(Umean, y) / H
    
    # Print Ubulk value
    print(f"Ubulk: {Ubulk}")
    
    # Plotting
    plt.figure(figsize=(12, 10))
    plt.plot(Umean, y, marker='o', linestyle='-', color='b', label=r"Moser1999, $Re_{\tau}=395$")
    
    # Labels
    plt.xlabel(r'$\overline{U}$', fontsize=14, fontname='serif')
    plt.ylabel(r'$y/H$', fontsize=14, fontname='serif')
    
    # Set y-axis range
    plt.ylim(0.0, 1.0)
    
    # Turn off grid
    plt.grid(False)
    
    # Legend
    plt.legend(fontsize=12)
    
    # Use LaTeX interpreter for the figure
    plt.rcParams.update({'mathtext.default': 'regular'})
    
    # Set font to serif for all text elements
    plt.rcParams['font.family'] = 'serif'
    
    # Set global font size
    plt.rcParams.update({'font.size': 14})
    
    # Show plot
    plt.tight_layout()
    plt.show()
    
except FileNotFoundError:
    print(f"The file at {file_path} does not exist.")
except pd.errors.EmptyDataError:
    print("The file is empty.")
except pd.errors.ParserError:
    print("Error parsing the file.")

