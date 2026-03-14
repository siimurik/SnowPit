import pandas as pd

def create_inverted_smet(csv_input, smet_output):
    # Load the CSV
    df = pd.read_csv(csv_input)

    # Unit Conversions
    df['TA_K'] = df['Temp_C'] + 273.15            # Real Air Temp
    df['TSG_K'] = df['Soil_Temp_320cm'] + 273.15    # Real Ground Temp
    df['TSS_K'] = int(-999)                          # Fixed 0°C Surface
    df['RH_frac'] = df['RH_%'] / 100.0
    
    smet_data = []

    # --- 1. Create the 2024-03-31T23:00 Anchor Line ---
    # This prevents SNOWPACK from crashing when it looks for data at the 23:50 start time
    first_row = df.iloc[0]
    anchor_line = (
        f"2024-03-31T23:00 "
        f"{first_row['TA_K']:.2f} "
        f"{first_row['RH_frac']:.3f} "
        f"{first_row['Air_Vel_m/s_10m']:.2f} "
        f"0.0 "
        f"0.00000 "
        f"{first_row['TSG_K']:.2f} "
#        f"{first_row['TSS_K']:.2f} "
        f"-999 "
        f"-999 "
    )
    smet_data.append(anchor_line)

    # --- 2. Process the rest of the data ---
    for _, row in df.iterrows():
        line = (
            f"{row['Time']} "            
            f"{row['TA_K']:.2f} "        
            f"{row['RH_frac']:.3f} "     
            f"{row['Air_Vel_m/s_10m']:.2f} " 
            f"0.0 "                      
            f"{row['Prec_m/h']:.5f} "    
            f"{row['TSG_K']:.2f} "       
            #f"{row['TSS_K']:.2f} "      
            f"-999 "
            f"-999"
        )
        smet_data.append(line)

    header = """SMET 1.1 ASCII
[HEADER]
station_id       = snow_storage
station_name     = Snow Storage Inverted
latitude         = 59.3980
longitude        = 24.6027
altitude         = 33.16
nodata           = -999
tz               = 0
fields           = timestamp TA RH VW ISWR PSUM TSG TSS HS
[DATA]
"""
    with open(smet_output, 'w') as f:
        f.write(header + "\n".join(smet_data) + "\n")
    
    print(f"Created {smet_output} with anchor line at 2024-03-31T23:00")

if __name__ == "__main__":
    create_inverted_smet('DATA_2024.csv', 'snow_storage.smet')