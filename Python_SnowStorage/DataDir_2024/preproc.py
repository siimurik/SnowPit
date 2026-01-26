import pandas as pd
from datetime import datetime
import numpy as np

# Load the Excel files - skip first 2 rows which contain headers/metadata
print("Loading data files...")
weather_df = pd.read_excel('Tallinn-Harku-2004-2024.xlsx', skiprows=2)
soil_df = pd.read_excel('Soil_temperatures_data_Tallinn.xlsx', skiprows=1)

# Assign proper column names based on the file structure
weather_df.columns = ['Year', 'Month', 'Day', 'Time', 'Glo_Sol_Ir_W/m2', 
                      'Pressure_Sea', 'Pressure_Station', 'Prec_mm/h', 'RH_%', 
                      'Temp_C', 'Temp_Min', 'Temp_Max', 'Wind_Dir', 
                      'Air_Vel_m/s_10m', 'Wind_Max']

soil_df.columns = ['Year', 'Month', 'Day', 'Time', 'Soil_20cm', 'Soil_40cm', 
                   'Soil_80cm', 'Soil_120cm', 'Soil_160cm', 'Soil_240cm', 'Soil_320cm']

print(f"Weather data loaded: {len(weather_df)} rows")
print(f"Soil data loaded: {len(soil_df)} rows")

# Select only needed columns
weather_filtered = weather_df[['Year', 'Month', 'Day', 'Time', 'Temp_C', 
                                'Air_Vel_m/s_10m', 'Prec_mm/h', 'Glo_Sol_Ir_W/m2', 'RH_%']].copy()
soil_filtered = soil_df[['Year', 'Month', 'Day', 'Time', 'Soil_320cm']].copy()
soil_filtered.rename(columns={'Soil_320cm': 'Soil_Temp_320cm'}, inplace=True)

print("\nProcessing datetime columns...")

# Create datetime column for both dataframes
def create_datetime(df):
    # Convert to numeric types
    df.loc[:, 'Year'] = pd.to_numeric(df['Year'], errors='coerce')
    df.loc[:, 'Month'] = pd.to_numeric(df['Month'], errors='coerce')
    df.loc[:, 'Day'] = pd.to_numeric(df['Day'], errors='coerce')
    
    # Drop rows with NaN in date columns
    df = df.dropna(subset=['Year', 'Month', 'Day']).copy()
    
    # Handle time format
    df.loc[:, 'Time'] = df['Time'].astype(str).str.strip()
    
    # Extract hour
    df.loc[:, 'Hour'] = df['Time'].str.split(':').str[0]
    df.loc[:, 'Hour'] = pd.to_numeric(df['Hour'], errors='coerce').fillna(0).astype(int)
    
    # Create datetime
    df.loc[:, 'DateTime'] = pd.to_datetime(
        df['Year'].astype(int).astype(str) + '-' +
        df['Month'].astype(int).astype(str).str.zfill(2) + '-' +
        df['Day'].astype(int).astype(str).str.zfill(2) + 'T' +
        df['Hour'].astype(str).str.zfill(2) + ':00',
        errors='coerce'
    )
    
    # Drop rows where datetime couldn't be created
    df = df.dropna(subset=['DateTime']).copy()
    
    return df

weather_filtered = create_datetime(weather_filtered)
soil_filtered = create_datetime(soil_filtered)

print(f"Weather data after cleaning: {len(weather_filtered)} rows")
print(f"Soil data after cleaning: {len(soil_filtered)} rows")

# Filter date range: 2024-04-01 to 2024-08-01
print("\nFiltering date range (2024-04-01 to 2024-08-01)...")
start_date = datetime(2024, 4, 1)
end_date = datetime(2024, 8, 1, 23, 59, 59)

weather_filtered = weather_filtered[
    (weather_filtered['DateTime'] >= start_date) & 
    (weather_filtered['DateTime'] <= end_date)
].copy()

soil_filtered = soil_filtered[
    (soil_filtered['DateTime'] >= start_date) & 
    (soil_filtered['DateTime'] <= end_date)
].copy()

print(f"Weather records in date range: {len(weather_filtered)}")
print(f"Soil records in date range: {len(soil_filtered)}")

# Merge weather and soil data
print("\nMerging datasets...")
merged_df = pd.merge(
    weather_filtered[['DateTime', 'Temp_C', 'Air_Vel_m/s_10m', 'Prec_mm/h', 
                      'Glo_Sol_Ir_W/m2', 'RH_%']],
    soil_filtered[['DateTime', 'Soil_Temp_320cm']],
    on='DateTime',
    how='left'
)

# Convert numeric columns
for col in ['Temp_C', 'Air_Vel_m/s_10m', 'Prec_mm/h', 'Glo_Sol_Ir_W/m2', 'RH_%', 'Soil_Temp_320cm']:
    merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce')

# Convert precipitation from mm/h to m/h
merged_df['Prec_m/h'] = merged_df['Prec_mm/h'] / 1000.0

# Format datetime as string
merged_df['Time'] = merged_df['DateTime'].dt.strftime('%Y-%m-%dT%H:%M')

# Select final columns in the correct order
final_df = merged_df[['Time', 'Temp_C', 'Air_Vel_m/s_10m', 'Prec_m/h', 
                       'Glo_Sol_Ir_W/m2', 'RH_%', 'Soil_Temp_320cm']].copy()

# Replace NaN values with 0.0 where appropriate
final_df = final_df.fillna(0.0)

# Save to CSV
output_file = 'DATA_2024.csv'
final_df.to_csv(output_file, index=False)

print("\n" + "="*80)
print(f"✓ Data successfully assembled and saved to {output_file}")
print(f"✓ Date range: {final_df['Time'].min()} to {final_df['Time'].max()}")
print(f"✓ Total records: {len(final_df)}")
print(f"\nStatistics:")
print(f"  Temperature range: {final_df['Temp_C'].min():.1f}°C to {final_df['Temp_C'].max():.1f}°C")
print(f"  Total precipitation: {final_df['Prec_m/h'].sum():.4f} m")
print(f"  Avg wind speed: {final_df['Air_Vel_m/s_10m'].mean():.2f} m/s")
print(f"\nFirst 10 rows:")
print(final_df.head(10))
print(f"\nLast 5 rows:")
print(final_df.tail(5))