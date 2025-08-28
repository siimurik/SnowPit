"""
Upload:
>   git add snow.py
>   git commit -m "Your commit message here"
>   git push origin main
Download:
>   git pull origin main
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

def find_first_empty_cell_index(data, column_index):
    """
    Find the index of the first row where the specified column is empty.
    
    Parameters:
    data (list): List of rows
    column_index (int): Index of column to check
    
    Returns:
    int: Index of first empty cell, or None if not found
    """
    for index, row in enumerate(data):
        try:
            if index < len(data) and column_index < len(row) and row[column_index].strip() == '':
                return index
        except (IndexError, AttributeError):
            continue
    return None

def extract_column(data, column_index):
    """
    Extract a specific column from the data.
    
    Parameters:
    data (list): List of rows
    column_index (int): Index of column to extract
    
    Returns:
    list: List of values from the specified column
    """
    return [row[column_index] for row in data if column_index < len(row)]

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
                print(f"{i}       {float(value):.6e}")
            except ValueError:
                print(f"{i}       {value}")
        print("       ...")
        for i in range(length-5, length):
            value = data[i]
            try:
                print(f"{i}    {float(value):.6e}")
            except ValueError:
                print(f"{i}    {value}")
    else:
        for index, value in enumerate(data):
            try:
                print(f"{index}       {float(value):.6e}")
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

def solve_tdma(a, b, c, d, n):
    """
    Solves a tridiagonal system of equations using the Thomas algorithm.

    Parameters:
        a (list): Lower diagonal of the matrix.
        b (list): Main diagonal of the matrix.
        c (list): Upper diagonal of the matrix.
        d (list): Right-hand side vector.
        n (int): Size of the system.

    Returns:
        list: Solution vector.
    """
    c_prime = [0.0] * n
    d_prime = [0.0] * n
    x = [0.0] * n

    # Forward sweep
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]

    for i in range(1, n):
        c_prime[i] = c[i] / (b[i] - a[i] * c_prime[i - 1])
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / (b[i] - a[i] * c_prime[i - 1])

    # Back substitution
    x[-1] = d_prime[-1]

    for i in range(n - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]

    return x

def transient1D(t_o, h_o, d_ins, lam_i, D, dx=0.005, dt=10.0, h_i=99.75):
    """
    Temperature distribution, transient 1D, BC 3d, implicit method.

    Parameters:
        t_o (list): Outer temperature for each hour.
        h_o (list): Outer heat transfer coefficient for each hour.
        d_ins (float): Insulation thickness, m.
        lam_i (float): Thermal conductivity.
        D (float): Thermal diffusivity.
        dx (float): Cell size, m.
        dt (float): Time interval, s.
        h_i (float): Inner heat transfer coefficient, W/m2K.

    Returns:
        list: 2D list of temperature distribution for each node and hour.
    """
    t_i = 0.0  # Inner temperature, C
    delta_db = d_ins  # Layer thickness, m
    k_db = lam_i  # Thermal conductivity
    a_db = D  # Thermal diffusivity
    n_el = int(delta_db / dx)  # Number of numerical elements
    nr_hour = len(t_o)  # Number of hours
    nodes = n_el + 1  # Number of nodes

    # Initialize temperature distribution
    T_n = [0.0] * nodes
    T_nh = [[0.0 for _ in range(nr_hour)] for _ in range(nodes)]  # 2D list for results

    # Initialize arrays for Tri-Diagonal Matrix Algorithm (TDMA)
    a = [0.0] * nodes
    b = [0.0] * nodes
    c = [0.0] * nodes
    d = [0.0] * nodes

    nh = int(3600 / dt)  # Number of time intervals in one hour

    # Main loop over hours
    for h in range(nr_hour):
        dFo = a_db * dt / dx**2
        dBio_i = h_i * dx / k_db
        dBio_o = h_o[h] * dx / k_db

        # Time steps within one hour
        for _ in range(nh):
            # Set up the tridiagonal system
            b[0] = (1.0 + 2.0 * dFo + 2.0 * dFo * dBio_o)
            c[0] = -2.0 * dFo
            d[0] = T_n[0] + 2.0 * dFo * dBio_o * t_o[h]

            for j in range(1, nodes - 1):
                a[j] = -dFo
                b[j] = (1.0 + 2.0 * dFo)
                c[j] = -dFo
                d[j] = T_n[j]

            a[-1] = -2.0 * dFo
            b[-1] = (1.0 + 2.0 * dFo + 2.0 * dFo * dBio_i)
            d[-1] = T_n[-1] + 2.0 * dFo * dBio_i * t_i

            # Solve the tridiagonal system
            T_n = solve_tdma(a, b, c, d, nodes)

        # Store the result for this hour
        for j in range(nodes):
            T_nh[j][h] = T_n[j]

    return T_nh

def export_large_matrix(vectors, filename, value_delimiter=' ', vector_delimiter='\n', buffer_size=1000):
    """
    Efficiently exports large vectors as a matrix to a file (PyPy compatible).
    
    Args:
        vectors (list of lists): Matrix where each inner list is a row vector
        filename (str): Output file path
        value_delimiter (str): Separator between values in a row (default: space)
        vector_delimiter (str): Separator between rows (default: newline)
        buffer_size (int): Number of rows to write at once (memory optimization)
        
    Returns:
        bool: True if successful, False if error occurred
    """
    try:
        # Check if vectors list is empty
        if not vectors:
            raise ValueError("Empty vectors list provided")
        
        # Get expected length from first vector
        expected_length = len(vectors[0])
        
        # Check all vectors have the same length
        for i, vector in enumerate(vectors):
            if len(vector) != expected_length:
                raise ValueError(f"All vectors must have same length. Vector at index {i} has length {len(vector)}, expected {expected_length}")
        
        with open(filename, 'w') as f:
            # Process vectors in chunks to avoid high memory usage
            for i in range(0, len(vectors), buffer_size):
                chunk = vectors[i:i+buffer_size]
                
                # Build chunk content
                chunk_content = []
                for vector in chunk:
                    # Efficient string joining for large vectors
                    row = value_delimiter.join(str(x) for x in vector)
                    chunk_content.append(row)
                
                # Write chunk with proper delimiters
                chunk_str = vector_delimiter.join(chunk_content)
                f.write(chunk_str)
                
                # Add delimiter between chunks if needed
                if i + buffer_size < len(vectors):
                    f.write(vector_delimiter)
        
        return True
    
    except (IOError, OSError) as e:
        print(f"File error: {e}")
        return False
    except ValueError as e:
        print(f"Data consistency error: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error: {e}")
        return False

from datetime import datetime

def convert_datetime_to_unix(time_strings):
    """Convert ISO format datetime strings to Unix timestamps with error handling."""
    unix_timestamps = []
    for i, time_str in enumerate(time_strings):
        try:
            dt = datetime.fromisoformat(time_str)
            unix_timestamp = dt.timestamp()
            unix_timestamps.append(unix_timestamp)
        except ValueError as e:
            print(f"Warning: Could not parse time string at index {i}: '{time_str}' - {e}")
            # Option: append NaN or None, or skip entirely
            unix_timestamps.append(float('nan'))
    return unix_timestamps

def main():
    # Read CSV data (keeping first 18 columns)
    data = read_csv_with_encoding('Snow_Storage_Data.csv', columns_to_keep=list(range(18)))
    
    # Display first few values from YEAR column (index 9)
    print("First 5 values from YEAR column:")
    year_values = extract_column(data, 9)[:5]
    for value in year_values:
        print(value)
    
    # Find first empty cell in YEAR.1 column (assuming it's also at index 9)
    first_nan_index = find_first_empty_cell_index(data, 9)
    print(f"\nFirst empty cell index in YEAR.1 column: {first_nan_index}")

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

    #printVec(time_column)
    rdata = data[period_start_in:period_end_in+1] # rdata - "ranged data"
    #printAny(rdata, column_name="All data")

    # Get air temperature data and assign it to a vector
    air_temp_vec_raw = [row[7] for row in rdata]
    air_temp_vec = convert_to_type(air_temp_vec_raw, dtype=float)
    #printVec(air_temp_vec, column_name="Air temperature (Celsius)")

    # Get air velocity data and assign it to a vector
    air_vel_vec_raw = [row[15] for row in rdata]
    air_vel_vec = convert_to_type(air_vel_vec_raw, dtype=float)
    #printVec(air_vel_vec, column_name="Air velocity (m/s)")

    # Extract the amount of precipitation column from the data
    prec_vec_raw = [row[5] for row in rdata]
    prec_vec = convert_to_type(prec_vec_raw, dtype=float)
    #printVec(prec_vec, column_name="Precipitation (m/h)")

    # Extract the global solar irradiance (W/m2) column from the data
    glob_solir_vec_raw = [row[12] for row in rdata]
    glob_solir_vec = convert_to_type(glob_solir_vec_raw, dtype=float)
    #printVec(glob_solir_vec, column_name="Global solar irradiance (W/m2)")

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
    #printVec(T_sol_air_vec, column_name="Solar-Air (Celsius)")

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
    #printVec(Q_surf_vec, column_name="Surface power (W)")

    #------------------------------------------------------------------------

    L_f = 333.4E03 # J/kg; latent heat of fusion
    rho_snow = 411.0 # kg/m^3; density of snow
    
    # The rate of melted snow from surface melt
    #f_srf_melt_vec = []
    #for i in range(len(Q_surf_vec)):
    #    f_srf_melt_vec.append(Q_surf_vec[i]/(L_f * rho_snow)) # m^3/s

    f_srf_melt_vec = [Q_surf / (L_f*rho_snow) for Q_surf in Q_surf_vec]
    #printVec(f_srf_melt_vec, column_name="Surface melt rate (m^3/s)") 

    #------------------------------------------------------------------------

    # Hourly rate of melted snow from surface melt
    #hrly_srf_total_vec = []
    #for i in range(len(f_srf_melt_vec)):
    #    hrly_srf_total_vec.append(f_srf_melt_vec[i] * 3600.0) # m^3/h

    hrly_srf_total_vec = [f_srf_melt*3600.0 for f_srf_melt in f_srf_melt_vec]
    #printVec(hrly_srf_total_vec, column_name="Hourly total SMR (m^3/h)")

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
    #printVec(q_rain_vec, column_name="Hourly q rain (W/m^2)")

    #------------------------------------------------------------------------

    # Hourly rain volume
    #v_rain_vec = []
    #for i in range(len(air_temp_vec)):
    #    v_rain_vec.append(prec_vec[i] * A_surf * rho_water * c_water * air_temp_vec[i] / (L_f * rho_snow)) # m^3/h
    
    v_rain_vec = [prec * A_surf * rho_water * c_water * air_temp / (L_f * rho_snow)
                  for prec, air_temp in zip(prec_vec, air_temp_vec)]
    #printVec(v_rain_vec, column_name="Hourly melt rate from rain (m^3/h)")

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
    #printVec(SMR_temp_vec, column_name="SMR due to T")

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
    #printVec(SMR_rain_vec, column_name="SMR due to rain (m^3/h)")

    #------------------------------------------------------------------------


    #SMR_total_vec = []
    #for i in range(len(SMR_temp_vec)):
    #    SMR_total_vec.append(SMR_temp_vec[i] + SMR_rain_vec[i])

    # Efficient version with list comprehension
    SMR_total_vec = [
        SMR_temp + SMR_rain for SMR_temp, SMR_rain 
        in zip(SMR_temp_vec, SMR_rain_vec) 
    ]
    #printVec(SMR_total_vec, column_name="Combined total SMR (m^3/h)")

    #------------------------------------------------------------------------

    # Cumulative sum of the SMR with rain together with temperature
    SMR_rainT_vec = cumsum(SMR_total_vec)
    #printVec(SMR_rainT_vec, column_name="Rain and T cumulative (m^3/h)")

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
    #printVec(emp1_SMR_vec, column_name="Empirical 1")

    #------------------------------------------------------------------------

    Psat_vec = []
    for i in range(len(air_temp_vec)):
        Psat_vec.append(Psat_WV(air_temp_vec[i] + 273.15)/10.0) # hPa; 100/1000 to convert to hPa
    #printVec(Psat_vec, column_name=" Water vapour saturation pressure (hPa)")

    # Extract the amount of RH column from the data
    RH_vec_raw = [row[17] for row in rdata]
    RH_vec = convert_to_type(RH_vec_raw, dtype=float)
    printVec(RH_vec, column_name="Relative Humidity (m/h)")
    
    # Checking to find the right range
    #printAny(time_column[period_start_in:period_end_in+1], column_name="Time")

    # Water steam pressure
    Pw_vec = []
    for i in range(len(RH_vec)):
        Pw_vec.append(Psat_vec[i]*RH_vec[i]/100.0) # kPa
    #printVec(Pw_vec, column_name = "Water steam pressure (kPa)")

    # Absolute humidity
    w_vec = []
    for i in range(len(Pw_vec)):
        w_vec.append(2.16679*Pw_vec[i]*1000/(273.15+air_temp_vec[i])) # kPa; 1000 to convert to kPa
    #printVec(w_vec, column_name="Absolute humidity (kPa)")

    # Surface melt rate from insulation thickness, air velocity, light intensity, air temperature and air humidity
    emp2_SMR_vec = []
    for i in range(len(air_temp_vec)):
        emp2_SMR_vec.append(- 0.97 - 0.097*(d_ins*100) + 0.164*air_vel_vec[i] + 0.00175*glob_solir_vec[i] 
                            + 0.102*air_temp_vec[i] + 0.192*w_vec[i]) # kg/m2/h
    #printVec(emp2_SMR_vec, column_name="Empirical 2")

    # Create a vector for emp1_SMR with condition ('wc')
    # This vector sets values to 0 if either the SMR or air temperature is below 0,
    # otherwise it retains the original SMR value from emp1_SMR_vec.
    emp1_SMR_wc_vec = [ # 'wc' - with condition
        0.0 if (smr < 0 or air_temp < 0) else smr 
        for smr, air_temp in zip(emp1_SMR_vec, air_temp_vec)
    ]
    #printVec(emp1_SMR_wc_vec, column_name="Empirical 1 (pos. cond)")

    # Create a vector for emp2_SMR_vec with condition ('cond')
    # This vector sets values to 0 if either the SMR or air temperature is below 0,
    # otherwise it retains the original SMR value from emp2_SMR_vec.
    emp2_SMR_wc_vec = [
        0.0 if (smr < 0 or air_temp < 0) else smr 
        for smr, air_temp in zip(emp2_SMR_vec, air_temp_vec)
    ]
    #printVec(emp2_SMR_wc_vec, column_name="Empirical 2 (pos. cond)")

    # Cumulative sum of pos. empirical methods
    emp1_SMR_wc_cs_vec = cumsum(emp1_SMR_wc_vec) # 'cs' - Cumulative Sum
    #printVec(emp1_SMR_wc_cs_vec, column_name="Cumulative sum - Emp 1")
    emp2_SMR_wc_cs_vec = cumsum(emp2_SMR_wc_vec) # 'cs' - Cumulative Sum
    #printVec(emp2_SMR_wc_cs_vec, column_name="Cumulative sum - Emp 2")

    # Reading in new data
    file_path = 'Tsi_Tso_Data.csv'
    data2 = read_csv_with_encoding(file_path)  # This will read all columns by default
    #printAny(data2)

    # If you want to keep specific columns
    #columns_to_keep = list(range(1))  # Define the indices of the columns you want to keep
    #data_with_specific_columns = read_csv_with_encoding(file_path, columns_to_keep)
    #printAny(data_with_specific_columns)

    Tsi_vec_raw = [row[0] for row in data2[1:]] # data2[1:] - exclude first row which contains the name
    #printAny(Tsi_vec_raw)
    Tsi_vec = convert_to_type(Tsi_vec_raw, dtype=float)
    #printVec(Tsi_vec, column_name="Tsi (Celsius)")

    Tso_vec_raw = [row[1] for row in data2[1:]] # data2[1:] - exclude first row which contains the name
    Tso_vec = convert_to_type(Tso_vec_raw, dtype=float)
    #printVec(Tso_vec, column_name="Tso (Celsius)")

    ho_vec = [
        6.0 + 4.0*vel if vel <= 5.0 else 7.41*(vel**0.78)
        for vel in air_vel_vec
    ]
    #printVec(ho_vec, column_name="# Air velocity (with cond)")

    # More constants for transient 1D solver
    h_i = 99.75  # W/m^2*K
    c_wet = 2.59E03 # J/(kg*K)
    rho_dry = 100.0 # kg/m^3
    moist_cont = 50.0 # %
    rho_wet = rho_dry + moist_cont/100.0*1000 # kg/m^3
    c_dry = 0.99E03 # J/(kg*K)
    c_water = 4.19E03 # J/(kg * K)

    c_wet = (1.0 - moist_cont/100.0)*c_dry + moist_cont/100.0*c_water
    #print(f"Specific heat capacity (wet): {c_wet/1000.0:.4} kJ/(kg*K).")

    D = lam_i/(c_wet * rho_wet) # m^2/s
    #print(f"Thermal diffusivity of the insulating material is {D:.4e} m^2/s.\n")

    # Snow outer and inner layer temperatures
    t_o = T_sol_air_vec
    h_o = ho_vec

    printVec(t_o)

    # Filtered time column of the observable period
    filt_Time = time_column[period_start_in:period_end_in+1]
    # Usage
    unix_timestamps = convert_datetime_to_unix(filt_Time)
    printVec(unix_timestamps, column_name="UNIX")
    print("")

    # Create matrix by pairing corresponding elements
    #scaled_RH = [rh * 0.01 for rh in RH_vec] # Example: RH goes from 72 to 0.72
    export_matrix = [[t, rh, to, ho] for t, rh, to, ho in zip(unix_timestamps, 
                            RH_vec, t_o, h_o)] # column format 

    # Export with tab separation (good for columnar data)
    export_large_matrix(export_matrix, "t_rh_to_ho.csv", value_delimiter=",")
    #export_large_matrix(export_matrix, "t_o_and_h_o.tsv", value_delimiter="\t")

    t_o_range = transient1D(t_o, h_o, d_ins, lam_i, D, dx=0.005, dt=10.0, h_i=h_i)
    #print(len(t_o_range[0]))

    # Print the temperature distribution for the first and last 5 hours
    for h in range(5):
        print(f"Hour {h}: Inner Temp = {t_o_range[-1][h]:.4e}, Outer Temp = {t_o_range[0][h]:.4e}")
    print("...")
    for h in range(-5, 0):
        print(f"Hour {len(ho_vec)-h-1}: Inner Temp = {t_o_range[-1][h]:.4e}, Outer Temp = {t_o_range[0][h]:.4e}")
    #printVec(t_o_range[0, :], name="Outer temperature of snow (°C)")  
    #printVec(t_o_range[-1, :], name="Internal temperature of snow (°C)")

    # Correction
    Tsi_vec = t_o_range[-1]
    Tso_vec = t_o_range[0]

    # Heat flux in
    qi_vec = []
    for i in range(len(Tsi_vec)):
        qi_vec.append((Tsi_vec[i] - 0.0)*h_i) # W/m^2
    #printVec(qi_vec, column_name="Heat flux in (W/m^2)")

    # Heat flux out
    qo_vec = []
    for i in range(len(T_sol_air_vec)):
        qo_vec.append((T_sol_air_vec[i] - Tso_vec[i]) * ho_vec[i]) # W/m^2
    #printVec(qo_vec, column_name="Heat flux out (W/m^2)")

    v_pc_vec = []
    for i in range(len(qi_vec)):
        v_pc_vec.append(qi_vec[i]/(L_f * rho_snow)) # m^3/(m^2*s)
    #printVec(v_pc_vec, column_name="Speed of phase change (m^3/(m^2*s))")

    v_pc_hourly_vec = []
    for i in range(len(v_pc_vec)):
        v_pc_hourly_vec.append(v_pc_vec[i] * 3600.0)  # m/h
    #printVec(v_pc_hourly_vec, column_name="Hourly speed of phase change (m^3/(m^2*h))")

    # Calculate hourly melt rate from solar heat flux
    #df_hfmr = pd.Series(np.where(df_air_temp > 0,  df_v_pc_hourly * rho_snow, 0.0))
    hfmr_vec = [
        v_pc_hourly * rho_snow if air_temp > 0 else 0.0
        for air_temp, v_pc_hourly in zip(air_temp_vec, v_pc_hourly_vec)
    ]
    #printVec(hfmr_vec, column_name="Hourly melt rate from solar heat flux")

    hfmr_cumsum_vec = cumsum(hfmr_vec)
    #printVec(hfmr_cumsum_vec, column_name="Cumulative hourly melt rate from solar heat flux")

    #rain_solar_hf_vec = []
    #for i in range(len(qo_vec)):
    #    rain_solar_hf_vec.append(q_rain_vec[i] + qo_vec[i]) # W/m^2
    rain_solar_hf_vec = [q_rain + qo for q_rain, qo in zip(q_rain_vec, qo_vec)]
    #printVec(rain_solar_hf_vec, column_name="Heat flux from rain and sun (W/m^2)")

    #wind_solar_rain_vec = []
    #for i in range(len(rain_solar_hf_vec)):
    #    wind_solar_rain_vec.append(rain_solar_hf_vec[i] / (T_sol_air_vec[i]-Tso_vec[i]))
    wind_solar_rain_vec = [rain_solar_hf/(T_sol_air-Tso) for rain_solar_hf, T_sol_air, Tso in zip(rain_solar_hf_vec, T_sol_air_vec, Tso_vec)]
    #printVec(wind_solar_rain_vec, column_name="Heat flux from wind, solar and rain (W/m^2)")

# End of main() section

if __name__ == "__main__":
    t0 = time.time()
    main()
    t1 = time.time()
    print(f"\nElapsed time: {t1-t0:.3e} sec.")
