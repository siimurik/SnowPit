import csv
import matplotlib.pyplot as plt

# Initialize lists to store the data
times = []
top_temps = []
bottom_temps = []

# Read the CSV file
with open('output_snow/temperature_results.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    
    current_time = None
    first_temp = None
    last_temp = None
    
    for row in csvreader:
        # Convert string values to appropriate data types
        time = float(row[0])
        node_index = int(row[1])
        position = float(row[2])
        temperature = float(row[3])
        
        # When we encounter a new timestep
        if current_time is None or time != current_time:
            # Save previous timestep's data if it exists
            if current_time is not None:
                times.append(current_time)
                top_temps.append(first_temp)
                bottom_temps.append(last_temp)
            
            # Start new timestep
            current_time = time
            first_temp = temperature  # First row is top layer
            last_temp = temperature   # Initialize last_temp, will be overwritten
        
        # Always update last_temp (will end up with last row's value)
        last_temp = temperature
    
    # Add the last timestep's data
    if current_time is not None:
        times.append(current_time)
        top_temps.append(first_temp)
        bottom_temps.append(last_temp)

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(times, top_temps, label='Top Layer Temperature (°C)')
plt.plot(times, bottom_temps, label='Bottom Layer Temperature (°C)')

# Add labels and title
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
plt.title('Temperature Evolution at Top and Bottom Layers')
plt.legend()
plt.grid(True)

# Show the plot
plt.tight_layout()
plt.show()