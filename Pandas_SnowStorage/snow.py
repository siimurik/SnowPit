"""
>   git add snow.py
>   git commit -m "Your commit message here"
>   git push origin main
"""

import chardet
import math
import time
import csv

def detect_encoding(file_path):
    """
    Detect the encoding of the given CSV file.

    Parameters:
    file_path (str): The path to the CSV file.

    Returns:
    str: The detected encoding of the file.
    """
    with open(file_path, 'rb') as rawdata:
        result = chardet.detect(rawdata.read(100000))
        return result['encoding']

def read_csv_with_encoding(file_path, columns_to_keep=None):
    """
    Read the CSV file with the detected encoding and keep specified columns.

    Parameters:
    file_path (str): The path to the CSV file.
    columns_to_keep (list, optional): A list of indices of columns to keep. 
                                      If None, all columns will be kept.

    Returns:
    list: A list of rows with only the specified columns (or all columns if none specified).
    """
    encoding = detect_encoding(file_path)
    print(f"\nDetected encoding: {encoding}")

    data = []
    with open(file_path, 'r', encoding=encoding) as file:
        reader = csv.reader(file)
        header = next(reader)  # Read the header to determine the number of columns
        max_columns = len(header)

        # If no columns specified, keep all columns
        if columns_to_keep is None:
            columns_to_keep = list(range(max_columns))

        for row in reader:
            selected_row = [row[i] for i in columns_to_keep]
            data.append(selected_row)
    
    return data

# Function to find the index of a specific value in a given column
def find_index(column_data, value):
    for index, item in enumerate(column_data):
        if item == value:
            return index
    return None

def convert_to_type(data, dtype=float):
    if dtype not in [float, int]:
        raise ValueError("dtype must be either 'float' or 'int'")
    
    converted_data = []
    for value in data:
        try:
            converted_data.append(dtype(value))
        except ValueError:
            print(f"Warning: Could not convert '{value}' to {dtype}")
            converted_data.append(None)  # Optionally, handle conversion errors
    
    return converted_data

def printAny(data, column_name="Unnamed"):
    print("")
    length = len(data)
    
    # Determine the data type
    data_type = type(data[0]).__name__ if length > 0 else "Unknown"

    if length > 10:
        for i in range(5):
            print(f"{i}       {data[i]}")
        print("       ...")
        for i in range(length-5, length):
            print(f"{i}    {data[i]}")
    else:
        for index, value in enumerate(data):
            print(f"{index}       {value}")
    
    print(f"Name: {column_name}, Length: {length}, dtype: {data_type}")


def printVec(data, column_name="Unnamed"):
    print("")  # Optional: print a newline for formatting
    length = len(data)
    
    # Determine the data type
    data_type = type(data[0]).__name__ if length > 0 else "Unknown"

    if length > 10:
        for i in range(5):
            value = data[i]
            try:
                print(f"{i}       {float(value):.6f}")
            except ValueError:
                print(f"{i}       {value}")
        print("       ...")
        for i in range(length-5, length):
            value = data[i]
            try:
                print(f"{i}    {float(value):.6f}")
            except ValueError:
                print(f"{i}    {value}")
    else:
        for index, value in enumerate(data):
            try:
                print(f"{index}       {float(value):.6f}")
            except ValueError:
                print(f"{index}       {value}")

    print(f"Name: {column_name}, Length: {length}, dtype: {data_type}")

def cumsum(values):
    cum_sum = []
    total = 0
    for value in values:
        total += value
        cum_sum.append(total)
    return cum_sum

def calculate_exp(x, terms=20):
    result = 1.0
    term = 1.0
    for n in range(1, terms):
        term *= x / n
        result += term
    return 

def Psat_WV(T_K):
    """
    Water vapour saturation pressure
    W. Wagner and A. Pruß:" The IAPWS Formulation 1995 for the
    Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use ",
    Journal of Physical and Chemical Reference Data, June 2002 ,Volume 31, Issue 2, pp.
    387535)

    Returns Saturation vapor pressure (hPa)
    """
    Tc = 647.096  # Critical temperature, K
    Pc = 220640   # Critical pressure, hPa
    
    C1 = -7.85951783
    C2 = 1.84408259
    C3 = -11.7866497
    C4 = 22.6807411
    C5 = -15.9618719
    C6 = 1.80122502
    
    teta = 1 - T_K / Tc
    
    x = Tc / T_K * (C1 * teta + C2 * teta ** 1.5 + C3 * teta ** 3 + C4 * teta ** 3.5 + C5 * teta ** 4 + C6 * teta ** 7.5)
    
    x = math.exp(x) * Pc
    
    return x

def main():
    # Detect the encoding of the file
    with open('Snow_Storage_Data.csv', 'rb') as rawdata:
        result = chardet.detect(rawdata.read(100000))
        encoding = result['encoding']

    print(encoding)

    # Create an empty list to store the data
    data = []

    # Define the indices of the columns you want to keep
    columns_to_keep = []
    for i in range(18):
        columns_to_keep.append(i)  

    # Read the CSV file with the detected encoding
    with open('Snow_Storage_Data.csv', 'r', encoding=encoding) as file:
        reader = csv.reader(file)
        header = next(reader)  # Skip the first row (header)
        for row in reader:
            selected_row = [row[i] for i in columns_to_keep]
            data.append(selected_row)

    # Display the first column of each row
    for row in data[:5]:
        print(row[9])   # row[9] - YEAR column

    # Find the index where the first NaN value appears in the 'YEAR.1' column (assuming it's the 5th column, change if needed)
    first_nan_index = None
    for index, row in enumerate(data):
        try:
            if row[9].strip() == '':
                first_nan_index = index
                break
        except IndexError:
            continue

    print(first_nan_index)

    # Get '2023-04-01T00:00' format time data from 7th column
    time_column = [row[6] for row in data]

    # Find the start and end indices
    period_start_in = find_index(time_column, '2023-04-01T00:00')
    period_end_in   = find_index(time_column, '2023-08-31T23:00')

    # Slice the data to include only the desired period
    #df_wo_nan_period = data[period_start_in:period_end_in + 1]  # +1 to include the end index

    # Display the sliced data
    #for row in data[period_start_in:period_end_in]:
    #    print(row[7])

    printVec(time_column)
    rdata = data[period_start_in:period_end_in+1] # rdata - "ranged data"
    printAny(rdata, column_name="All data")

    # Get air temperature data and assign it to a vector
    air_temp_vec_raw = [row[7] for row in rdata]
    air_temp_vec = convert_to_type(air_temp_vec_raw, dtype=float)
    printVec(air_temp_vec, column_name="Air temperature (Celsius)")

    # Get air velocity data and assign it to a vector
    air_vel_vec_raw = [row[15] for row in rdata]
    air_vel_vec = convert_to_type(air_vel_vec_raw, dtype=float)
    printVec(air_vel_vec, column_name="Air velocity (m/s)")

    # Extract the amount of precipitation column from the data
    prec_vec_raw = [row[5] for row in rdata]
    prec_vec = convert_to_type(prec_vec_raw, dtype=float)
    printVec(prec_vec, column_name="Precipitation (m/h)")

    # Extract the global solar irradiance (W/m2) column from the data
    glob_solir_vec_raw = [row[12] for row in rdata]
    glob_solir_vec = convert_to_type(glob_solir_vec_raw, dtype=float)
    printVec(glob_solir_vec, column_name="Global solar irradiance (W/m2)")

    #------------------------------------------------------------------------

    # Heat transfer coefficient at the external surface
    h = 22.7 # W/(m^2K)
    # Solar light absorptivity
    alpha = 0.8
    # Correction factor for horizontal surface
    T_cor_fact = 4.0 # °C
    
    #T_sol_air_vec = []
    #for i in range(len(air_temp_vec)):
    #    T_sol_air_vec.append(alpha * glob_solir_vec[i] / h + air_temp_vec[i] - T_cor_fact)
    
    # Alternative version:
    # Improved efficiency with list comprehension
    T_sol_air_vec = [
        alpha * glob_solir / h + air_temp - T_cor_fact 
        for glob_solir, air_temp in 
        zip(glob_solir_vec, air_temp_vec)
    ]
    printVec(T_sol_air_vec, column_name="Solar-Air (Celsius)")

    #------------------------------------------------------------------------

    # Insulation layer thickness
    d_ins = 0.1 # m
    # Thermal conductivity for the insulating material
    lam_i = 0.32 # W/(mK)
    # The surface area (m2) of the pile of snow
    A_surf = 210.0 # m^2
    
    # The rate of heat transfer from the snow pile to the air
    #Q_surf_vec = []
    #for i in range(len(T_sol_air_vec)):
    #    Q_surf_vec.append(A_surf * lam_i / d_ins * (T_sol_air_vec[i] - 0.0))    # W

    # More efficient version
    Q_surf_vec = [
        A_surf * lam_i / d_ins * (T_sol_air - 0.0)
        for T_sol_air in T_sol_air_vec
    ]
    printVec(Q_surf_vec, column_name="Surface power (W)")

    #------------------------------------------------------------------------

    L_f = 333.4E03 # J/kg; latent heat of fusion
    rho_snow = 411.0 # kg/m^3; density of snow
    
    # The rate of melted snow from surface melt
    #f_srf_melt_vec = []
    #for i in range(len(Q_surf_vec)):
    #    f_srf_melt_vec.append(Q_surf_vec[i]/(L_f * rho_snow)) # m^3/s

    f_srf_melt_vec = [Q_surf / (L_f*rho_snow) for Q_surf in Q_surf_vec]
    printVec(f_srf_melt_vec, column_name="Surface melt rate (m^3/s)") 

    #------------------------------------------------------------------------

    # Hourly rate of melted snow from surface melt
    #hrly_srf_total_vec = []
    #for i in range(len(f_srf_melt_vec)):
    #    hrly_srf_total_vec.append(f_srf_melt_vec[i] * 3600.0) # m^3/h

    hrly_srf_total_vec = [f_srf_melt*3600.0 for f_srf_melt in f_srf_melt_vec]
    printVec(hrly_srf_total_vec, column_name="Hourly total SMR (m^3/h)")

    #------------------------------------------------------------------------

    # Initialize q_rain_vec with zeros
    #q_rain_vec = [0.0] * len(air_vel_vec)

    # Constants
    rho_water = 1000.0  # kg/m3
    c_water   = 4.19E03  # J/(kg*K)

    # Initialize pos_temp_mark and calculate heat flux
    #pos_temp_mark = []
    #for i in range(len(air_temp_vec)):
    #    if air_temp_vec[i] > 0.0:
    #        q_rain_vec[i] = prec_vec[i] * rho_water * c_water * air_temp_vec[i] / 3600.0
    #        pos_temp_mark.append(True)
    #    else:
    #        pos_temp_mark.append(False)

    # Update q_rain_vec where air_temp is greater than 0.0
    #q_rain_vec = [
    #    prec * rho_water * c_water * air_temp / 3600.0
    #    for prec, air_temp in zip(prec_vec, air_temp_vec) 
    #    if air_temp > 0.0 
    #]

    # Initialize q_rain_vec with zeros
    q_rain_vec = [
        prec * rho_water * c_water * air_temp / 3600.0 if air_temp > 0.0 else 0.0
        for prec, air_temp in zip(prec_vec, air_temp_vec)
    ]

    # Display results
    #print(pos_temp_mark[:5])
    printVec(q_rain_vec, column_name="Hourly q rain (m^3/h)")

    #------------------------------------------------------------------------

    # Hourly rain volume
    #v_rain_vec = []
    #for i in range(len(air_temp_vec)):
    #    v_rain_vec.append(prec_vec[i] * A_surf * rho_water * c_water * air_temp_vec[i] / (L_f * rho_snow)) # m^3/h
    
    v_rain_vec = [prec * A_surf * rho_water * c_water * air_temp / (L_f * rho_snow)
                  for prec, air_temp in zip(prec_vec, air_temp_vec)]
    printVec(v_rain_vec, column_name="Hourly melt rate from rain (m^3/h)")

    #------------------------------------------------------------------------

    # Calculate the surface melt rate where the temperature is greater than 0
    #SMR_temp_vec = [0.0] * len(air_vel_vec)
    #for i in range(len(air_temp_vec)):
    #    if air_temp_vec[i] > 0.0:
    #        SMR_temp_vec[i] = hrly_srf_total_vec[i] * rho_snow / A_surf # m^3/h
    
    # Efficient version with list comprehension
    SMR_temp_vec = [
        hrly_srf_total*rho_snow/A_surf if air_temp > 0.0 else 0.0
        for hrly_srf_total, air_temp in zip(hrly_srf_total_vec, air_temp_vec)
    ]
    printVec(SMR_temp_vec, column_name="SMR due to T")

    #------------------------------------------------------------------------

    # Surface melt rate due to rain
    #SMR_rain_vec = []
    #for i in range(len(v_rain_vec)):
    #    SMR_rain_vec.append(v_rain_vec[i] * rho_snow / A_surf) # m^3/h

    # Efficient version with list comprehension
    SMR_rain_vec = [
        v_rain * rho_snow / A_surf
        for v_rain in v_rain_vec
    ]
    printVec(SMR_rain_vec, column_name="SMR due to rain (m^3/h)")

    #------------------------------------------------------------------------


    #SMR_total_vec = []
    #for i in range(len(SMR_temp_vec)):
    #    SMR_total_vec.append(SMR_temp_vec[i] + SMR_rain_vec[i])

    # Efficient version with list comprehension
    SMR_total_vec = [
        SMR_temp + SMR_rain for SMR_temp, SMR_rain 
        in zip(SMR_temp_vec, SMR_rain_vec) 
    ]
    printVec(SMR_total_vec, column_name="Combined toal SMR (m^3/h)")

    #------------------------------------------------------------------------

    # Cumulative sum of the SMR with rain together with temperature
    SMR_rainT_vec = cumsum(SMR_total_vec)
    printVec(SMR_rainT_vec, column_name="Rain and T cumulative (m^3/h)")

    #------------------------------------------------------------------------

    #emp1_SMR_vec = []
    #for i in range(len(air_temp_vec)):
    #    emp1_SMR_vec.append(-0.09 + 0.00014*glob_solir_vec[i] + 0.0575*air_temp_vec[i] + 
    #                        0.0012*air_temp_vec[i]*air_vel_vec[i] - 0.18*air_temp_vec[i]*d_ins) # kg/m2/h

    # Efficient version with list comprehension
    emp1_SMR_vec = [
        -0.09 + 0.00014*glob_solir + 0.0575*air_temp + 0.0012*air_temp*air_vel - 0.18*air_temp*d_ins
        for glob_solir, air_temp, air_vel in zip(glob_solir_vec, air_temp_vec, air_vel_vec)
    ]
    printVec(emp1_SMR_vec, column_name="Empirical 1")

    #------------------------------------------------------------------------

    Psat_vec = []
    for i in range(len(air_temp_vec)):
        Psat_vec.append(Psat_WV(air_temp_vec[i] + 273.15)/10.0) # hPa; 100/1000 to convert to hPa
    printVec(Psat_vec, column_name=" Water vapour saturation pressure (hPa)")

    # Extract the amount of RH precipitation column from the data
    RH_perc_vec_raw = [row[17] for row in rdata]
    RH_perc_vec = convert_to_type(RH_perc_vec_raw, dtype=float)
    printVec(prec_vec, column_name="Relative Humidity Percipitation (m/h)")

    # Water steam pressure
    Pw_vec = []
    for i in range(len(RH_perc_vec)):
        Pw_vec.append(Psat_vec[i]*RH_perc_vec[i]/100.0) # kPa
    printVec(Pw_vec, column_name = "Water steam pressure (kPa)")

    # Absolute humidity
    w_vec = []
    for i in range(len(Pw_vec)):
        w_vec.append(2.16679*Pw_vec[i]*1000/(273.15+air_temp_vec[i])) # kPa; 1000 to convert to kPa
    printVec(w_vec, column_name="Absolute humidity (kPa)")

    # Surface melt rate from insulation thickness, air velocity, light intensity, air temperature and air humidity
    emp2_SMR_vec = []
    for i in range(len(air_temp_vec)):
        emp2_SMR_vec.append(- 0.97 - 0.097*(d_ins*100) + 0.164*air_vel_vec[i] + 0.00175*glob_solir_vec[i] 
                            + 0.102*air_temp_vec[i] + 0.192*w_vec[i]) # kg/m2/h
    printVec(emp2_SMR_vec, column_name="Empirical 2")

    # Create a vector for emp1_SMR with condition ('wc')
    # This vector sets values to 0 if either the SMR or air temperature is below 0,
    # otherwise it retains the original SMR value from emp1_SMR_vec.
    emp1_SMR_wc_vec = [ # 'wc' - with condition
        0.0 if (smr < 0 or air_temp < 0) else smr 
        for smr, air_temp in zip(emp1_SMR_vec, air_temp_vec)
    ]
    printVec(emp1_SMR_wc_vec, column_name="Empirical 1 (pos. cond)")

    # Create a vector for emp2_SMR_vec with condition ('cond')
    # This vector sets values to 0 if either the SMR or air temperature is below 0,
    # otherwise it retains the original SMR value from emp2_SMR_vec.
    emp2_SMR_wc_vec = [
        0.0 if (smr < 0 or air_temp < 0) else smr 
        for smr, air_temp in zip(emp2_SMR_vec, air_temp_vec)
    ]
    printVec(emp2_SMR_wc_vec, column_name="Empirical 2 (pos. cond)")

    # Cumulative sum of pos. empirical methods
    emp1_SMR_wc_cs_vec = cumsum(emp1_SMR_wc_vec) # 'cs' - Cumulative Sum
    printVec(emp1_SMR_wc_cs_vec, column_name="Cumulative sum - Emp 1")
    emp2_SMR_wc_cs_vec = cumsum(emp2_SMR_wc_vec) # 'cs' - Cumulative Sum
    printVec(emp2_SMR_wc_cs_vec, column_name="Cumulative sum - Emp 2")

    # Reading in new data
    file_path = 'Tsi_Tso_Data.csv'
    data2 = read_csv_with_encoding(file_path)  # This will read all columns by default
    printAny(data2)

    # If you want to keep specific columns
    #columns_to_keep = list(range(1))  # Define the indices of the columns you want to keep
    #data_with_specific_columns = read_csv_with_encoding(file_path, columns_to_keep)
    #printAny(data_with_specific_columns)

    Tsi_vec_raw = [row[0] for row in data2[1:]] # data2[1:] - exclude first row which contains the name
    #printAny(Tsi_vec_raw)
    Tsi_vec = convert_to_type(Tsi_vec_raw, dtype=float)
    printVec(Tsi_vec, column_name="Tsi (Celsius)")

    Tso_vec_raw = [row[1] for row in data2[1:]] # data2[1:] - exclude first row which contains the name
    Tso_vec = convert_to_type(Tso_vec_raw, dtype=float)
    printVec(Tso_vec, column_name="Tso (Celsius)")

    ho_vec = [
        6.0 + 4.0*vel if vel <= 5.0 else 7.41*(vel**0.78)
        for vel in air_vel_vec
    ]
    printVec(ho_vec, column_name="# Air velocity (with cond)")

    # Heat transfer coefficient at the internal surface:
    h_i = 99.75 # W/m^2*K
    # Heat flux in
    qi_vec = []
    for i in range(len(Tsi_vec)):
        qi_vec.append((Tsi_vec[i] - 0.0)*h_i) # W/m^2
    printVec(qi_vec, column_name="Heat flux in (W/m^2)")

    # Heat flux out
    qo_vec = []
    for i in range(len(T_sol_air_vec)):
        qo_vec.append((T_sol_air_vec[i] - Tso_vec[i]) * ho_vec[i]) # W/m^2
    printVec(qo_vec, column_name="Heat flux put (W/m^2)")

    v_pc_vec = []
    for i in range(len(qi_vec)):
        v_pc_vec.append(qi_vec[i]/(L_f * rho_snow)) # m^3/(m^2*s)
    printVec(v_pc_vec, column_name="Speed of phase change (m^3/(m^2*s))")

    v_pc_hourly_vec = []
    for i in range(len(v_pc_vec)):
        v_pc_hourly_vec.append(v_pc_vec[i] * 3600.0)  # m/h
    printVec(v_pc_hourly_vec, column_name="Hourly speed of phase change (m^3/(m^2*h))")

    # Calculate hourly melt rate from solar heat flux
    #df_hfmr = pd.Series(np.where(df_air_temp > 0,  df_v_pc_hourly * rho_snow, 0.0))
    hfmr_vec = [
        v_pc_hourly * rho_snow if air_temp > 0 else 0.0
        for air_temp, v_pc_hourly in zip(air_temp_vec, v_pc_hourly_vec)
    ]
    printVec(hfmr_vec, column_name="Hourly melt rate from solar heat flux")

    hfmr_cumsum_vec = cumsum(hfmr_vec)
    printVec(hfmr_cumsum_vec, column_name="Cumulative hourly melt rate from solar heat flux")

    #rain_solar_hf_vec = []
    #for i in range(len(qo_vec)):
    #    rain_solar_hf_vec.append(q_rain_vec[i] + qo_vec[i]) # W/m^2
    rain_solar_hf_vec = [q_rain + qo for q_rain, qo in zip(q_rain_vec, qo_vec)]
    printVec(rain_solar_hf_vec, column_name="Heat flux from rain and sun (W/m^2)")

    #wind_solar_rain_vec = []
    #for i in range(len(rain_solar_hf_vec)):
    #    wind_solar_rain_vec.append(rain_solar_hf_vec[i] / (T_sol_air_vec[i]-Tso_vec[i]))
    wind_solar_rain_vec = [rain_solar_hf/(T_sol_air-Tso) for rain_solar_hf, T_sol_air, Tso in zip(rain_solar_hf_vec, T_sol_air_vec, Tso_vec)]
    printVec(wind_solar_rain_vec, column_name="Heat flux from wind, solar and rain (W/m^2)")

# End of main() section

if __name__ == "__main__":
    t0 = time.time()
    main()
    t1 = time.time()
    print(f"\nElapsed time: {t1-t0:.3e} sec.")