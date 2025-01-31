import chardet
import math
import csv

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

# Function to find the index of a specific value in a given column
def find_index(column_data, value):
    for index, item in enumerate(column_data):
        if item == value:
            return index
    return 

# Find the start and end indices
period_start_in = find_index(time_column, '2023-04-01T00:00')
period_end_in   = find_index(time_column, '2023-08-31T23:00')

# Slice the data to include only the desired period
#df_wo_nan_period = data[period_start_in:period_end_in + 1]  # +1 to include the end index

# Display the sliced data
#for row in data[period_start_in:period_end_in]:
#    print(row[7])

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

# Heat transfer coefficient at the external surface
h = 22.7 # W/(m^2K)
# Solar light absorptivity
alpha = 0.8
# Correction factor for horizontal surface
T_cor_fact = 4.0 # °C
T_sol_air_vec = []
for i in range(len(air_temp_vec)):
    T_sol_air_vec.append(alpha * glob_solir_vec[i] / h + air_temp_vec[i] - T_cor_fact)
printVec(T_sol_air_vec, column_name="Solar-Air (Celsius)")

# Insulation layer thickness
d_ins = 0.1 # m
# Thermal conductivity for the insulating material
lam_i = 0.32 # W/(mK)
# The surface area (m2) of the pile of snow
A_surf = 210.0 # m^2
# The rate of heat transfer from the snow pile to the air
Q_surf_vec = []
for i in range(len(T_sol_air_vec)):
    Q_surf_vec.append(A_surf * lam_i / d_ins * (T_sol_air_vec[i] - 0.0))    # W
printVec(Q_surf_vec, column_name="Surface power (W)")

L_f = 333.4E03 # J/kg; latent heat of fusion
rho_snow = 411.0 # kg/m^3; density of snow
# The rate of melted snow from surface melt
f_srf_melt_vec = []
for i in range(len(Q_surf_vec)):
    f_srf_melt_vec.append(Q_surf_vec[i]/(L_f * rho_snow)) # m^3/s
printVec(f_srf_melt_vec, column_name="Surface melt rate (m^3/s)") 

# Hourly rate of melted snow from surface melt
hrly_srf_total_vec = []
for i in range(len(f_srf_melt_vec)):
    hrly_srf_total_vec.append(f_srf_melt_vec[i] * 3600) # m^3/h
printVec(hrly_srf_total_vec, column_name="Hourly rate of surface melt (m^3/h)")

# Initialize q_rain_vec with zeros
q_rain_vec = [0.0] * len(air_vel_vec)

# Constants
rho_water = 1000.0  # kg/m3
c_water   = 4.19E03  # J/(kg*K)

# Initialize pos_temp_mark and calculate heat flux
pos_temp_mark = []
for i in range(len(air_temp_vec)):
    if air_temp_vec[i] > 0.0:
        q_rain_vec[i] = prec_vec[i] * rho_water * c_water * air_temp_vec[i] / 3600.0
#        pos_temp_mark.append(True)
#    else:
#        pos_temp_mark.append(False)

# Display results
#print(pos_temp_mark[:5])
printVec(q_rain_vec, column_name="Hourly melt rate from surface (m^3/h)")

# Hourly rain volume
v_rain_vec = []
for i in range(len(air_temp_vec)):
    v_rain_vec.append(prec_vec[i] * A_surf * rho_water * c_water * air_temp_vec[i] / (L_f * rho_snow)) # m^3/h
printVec(v_rain_vec, column_name="Hourly melt rate from rain (m^3/h)")

# Calculate the surface melt rate where the temperature is greater than 0
SMR_temp_vec = [0.0] * len(air_vel_vec)
for i in range(len(air_temp_vec)):
    if air_temp_vec[i] > 0.0:
        SMR_temp_vec[i] = hrly_srf_total_vec[i] * rho_snow / A_surf # m^3/h
printVec(SMR_temp_vec, column_name="SMR due to T")

# Surface melt rate due to rain
SMR_rain_vec = []
for i in range(len(v_rain_vec)):
    SMR_rain_vec.append(v_rain_vec[i] * rho_snow / A_surf) # m^3/h
printVec(SMR_rain_vec, column_name="SMR due to rain (m^3/h)")

SMR_total_vec = []
for i in range(len(SMR_temp_vec)):
    SMR_total_vec.append(SMR_temp_vec[i] + SMR_rain_vec[i])
printVec(SMR_total_vec, column_name="Combined toal SMR (m^3/h)")

def cumsum(values):
    cum_sum = []
    total = 0
    for value in values:
        total += value
        cum_sum.append(total)
    return cum_sum

# Cumulative sum of the SMR with rain together with temperature
SMR_rainT_vec = cumsum(SMR_total_vec)
printVec(SMR_rainT_vec, column_name="Rain and T cumulative (m^3/h)")

emp1_SMR_vec = []
for i in range(len(air_temp_vec)):
    emp1_SMR_vec.append(-0.09 + 0.00014*glob_solir_vec[i] + 0.0575*air_temp_vec[i] + 
                        0.0012*air_temp_vec[i]*air_vel_vec[i] - 0.18*air_temp_vec[i]*d_ins) # kg/m2/h
printVec(emp1_SMR_vec, column_name="Empirical 1")

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
Psat_vec = []
for i in range(len(air_temp_vec)):
    Psat_vec.append(Psat_WV(air_temp_vec[i] + 273.15)/10.0) # hPa; 100/1000 to convert to hPa
printVec(Psat_vec, column_name=" Water vapour saturation pressure (hPa)")