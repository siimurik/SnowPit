#!/usr/bin/env python3
"""
Temperature Profile Comparison: FEM vs FDM
Compares temperature profiles from FEM and Implicit FDM simulations
"""

import csv
import matplotlib.pyplot as plt

# Configuration
FEM_DATA_PATH = 'output_snow/temperature_results.csv'
FDM_DATA_PATH = 'output_snow/combined_temp_data.csv'
PLOT_FIRST_N_POINTS = 100

# Plot styling
plt.style.use('seaborn-v0_8-darkgrid')
COLORS = {
    'fem_outer': '#1f77b4',      # Blue
    'fem_internal': '#d62728',   # Red
    'fdm_outer': '#17becf',      # Teal
    'fdm_internal': '#e377c2'    # Pink
}
LINE_STYLES = {
    'fem': '-',      # Solid line
    'fdm': '--'      # Dashed line
}
LINE_WIDTH = {
    'fem': 2.5,      # FEM line thickness
    'fdm': 2.0       # FDM line thickness
}

def read_fem_data(filepath):
    """Read FEM simulation data and extract outer/internal temperatures per time step."""
    times, outer_temps, internal_temps = [], [], []
    
    current_time = None
    first_temp = None
    last_temp = None
    
    with open(filepath, 'r') as f:
        for row in csv.reader(f):
            time = float(row[0])
            temperature = float(row[3])
            
            # New time step detected
            if current_time is None or time != current_time:
                # Save previous time step data
                if current_time is not None:
                    times.append(current_time)
                    outer_temps.append(first_temp)
                    internal_temps.append(last_temp)
                
                # Initialize new time step
                current_time = time
                first_temp = temperature
            
            # Update last temperature for current time step
            last_temp = temperature
    
    # Add final time step
    if current_time is not None:
        times.append(current_time)
        outer_temps.append(first_temp)
        internal_temps.append(last_temp)
    
    return times, outer_temps, internal_temps

def read_fdm_data(filepath):
    """Read FDM simulation data from CSV."""
    times, outer_temps, internal_temps = [], [], []
    
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header row
        
        for row in reader:
            times.append(float(row[0]))           # Time
            outer_temps.append(float(row[1]))     # Outer surface temp
            internal_temps.append(float(row[2]))  # Internal temp
    
    return times, outer_temps, internal_temps

def create_comparison_plot(fem_data, fdm_data, n_points=None):
    """Create comparison plot of FEM vs FDM temperature data."""
    # Unpack data
    fem_times, fem_outer, fem_internal = fem_data
    fdm_times, fdm_outer, fdm_internal = fdm_data
    
    # Use FDM time series for both datasets (this is the key!)
    # Limit FDM data to n_points first
    if n_points:
        plot_times = fdm_times[:n_points]
        fdm_outer = fdm_outer[:n_points]
        fdm_internal = fdm_internal[:n_points]
    else:
        plot_times = fdm_times
    
    # Limit FEM data to same number of points as FDM
    fem_outer = fem_outer[:len(plot_times)]
    fem_internal = fem_internal[:len(plot_times)]
    
    # Create plot
    plt.figure(figsize=(12, 7))
    
    # Plot BOTH datasets against the SAME time series (FDM times)
    plt.plot(plot_times, fem_outer, 
             color=COLORS['fem_outer'], linestyle=LINE_STYLES['fem'], 
             linewidth=LINE_WIDTH['fem'], label='Outer Layer (FEM)')
    plt.plot(plot_times, fem_internal, 
             color=COLORS['fem_internal'], linestyle=LINE_STYLES['fem'], 
             linewidth=LINE_WIDTH['fem'], label='Internal Layer (FEM)')
    
    # Plot FDM data - using same time series
    plt.plot(plot_times, fdm_outer, 
             color=COLORS['fdm_outer'], linestyle=LINE_STYLES['fdm'], 
             linewidth=LINE_WIDTH['fdm'], label='Outer Layer (Implicit FDM)')
    plt.plot(plot_times, fdm_internal, 
             color=COLORS['fdm_internal'], linestyle=LINE_STYLES['fdm'], 
             linewidth=LINE_WIDTH['fdm'], label='Internal Layer (Implicit FDM)')
    
    # Formatting
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Temperature (°C)', fontsize=12)
    plt.title('Temperature Profile Comparison: FEM vs Implicit FDM', fontsize=14)
    plt.legend(fontsize=10, loc='lower right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return plt.gcf()

# Main execution
def main():
    print("Reading temperature data...")
    
    # Read data from both sources
    fem_data = read_fem_data(FEM_DATA_PATH)
    fdm_data = read_fdm_data(FDM_DATA_PATH)
    
    print(f"FEM data points: {len(fem_data[0])}")
    print(f"FDM data points: {len(fdm_data[0])}")
    
    # Create and display plot
    fig = create_comparison_plot(fem_data, fdm_data, PLOT_FIRST_N_POINTS)
    plt.show()
    
    # Uncomment to save plot
    fig.savefig('fem_vs_fdm_comparison.png', dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    main()