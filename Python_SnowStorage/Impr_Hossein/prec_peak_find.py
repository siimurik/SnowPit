import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Read the CSV file
df = pd.read_csv('DATA_2024.csv')

# Convert Time to datetime
df['Time'] = pd.to_datetime(df['Time'])

# Convert precipitation from m/h to mm/h (multiply by 1000)
df['Prec_mm/h'] = df['Prec_m/h'] * 1000

# Set tolerance threshold
tolerance = 1.30  # mm/h

# Filter for values above tolerance
peaks = df[df['Prec_mm/h'] > tolerance].copy()

# Sort by precipitation in descending order to show highest peaks first
peaks_sorted = peaks.sort_values('Prec_mm/h', ascending=False)

# Display results
print(f"Precipitation Peak Analysis")
print(f"Tolerance threshold: {tolerance} mm/h")
print(f"\nNumber of peaks found: {len(peaks)}")

if len(peaks) > 0:
    print(f"\nHighest peaks above {tolerance} mm/h:\n")
    print(peaks_sorted[['Time', 'Prec_m/h', 'Prec_mm/h']].to_string(index=False))
    
    print(f"\nHighest precipitation value: {peaks_sorted['Prec_mm/h'].max():.2f} mm/h")
    print(f"Occurred at: {peaks_sorted.iloc[0]['Time']}")
else:
    print(f"\nNo precipitation values found above {tolerance} mm/h")

# Create the plot
fig, ax = plt.subplots(figsize=(12, 6))

# Plot all precipitation data
ax.plot(df['Time'], df['Prec_mm/h'], color='steelblue', linewidth=1.5, label='Precipitation')

# Highlight peaks above tolerance
if len(peaks) > 0:
    ax.scatter(peaks['Time'], peaks['Prec_mm/h'], 
              color='red', s=100, zorder=5, label=f'Peaks > {tolerance} mm/h', 
              edgecolors='darkred', linewidths=2)

# Add horizontal line for tolerance threshold
ax.axhline(y=tolerance, color='orange', linestyle='--', linewidth=2, 
          label=f'Tolerance ({tolerance} mm/h)')

# Format the plot
ax.set_xlabel('Time', fontsize=12)
ax.set_ylabel('Precipitation (mm/h)', fontsize=12)
ax.set_title('Precipitation Over Time with Peak Detection', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(loc='upper right')

# Format x-axis dates
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
plt.xticks(rotation=45, ha='right')

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Save and show the plot
plt.savefig('precipitation_peaks.png', dpi=300, bbox_inches='tight')
print("\nPlot saved as 'precipitation_peaks.png'")
plt.show()