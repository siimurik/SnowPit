import numpy as np
import pyarrow.parquet as pq

def printNP(arr, name=None):
    """
    Prints a NumPy array in the usual format with optional name, length, and dtype.
    """
    array_str = np.array2string(arr, separator=' ', threshold=6)
    length = arr.shape[0]
    dtype = arr.dtype
    if name:
        print(f"{array_str}\nName: {name}, Length: {length}, dtype: {dtype}")
    else:
        print(f"{array_str}\nLength: {length}, dtype: {dtype}")

def Psat_WV(T_K):
    """
    Water vapour saturation pressure.
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
    
    x = np.exp(x) * Pc
    
    return x

# Read the Parquet file using PyArrow
parquet_file = 'SnowStorageDATA.parquet'
table = pq.read_table(parquet_file)
data = table.to_pandas().to_numpy()

# Process the time column
time_column = data[:, 6].astype('datetime64[s]')
period_start_in = np.where(time_column == np.datetime64('2023-04-01T00:00'))[0][0]
period_end_in = np.where(time_column == np.datetime64('2023-08-31T23:00'))[0][0]
rdata = data[period_start_in:period_end_in + 1]

# Convert columns to NumPy arrays
air_temp_vec = rdata[:, 7].astype(float)
air_vel_vec = rdata[:, 15].astype(float)
prec_vec = rdata[:, 5].astype(float)
glob_solir_vec = rdata[:, 12].astype(float)

# Constants
h = 22.7
alpha = 0.8
T_cor_fact = 4.0
L_f = 333.4E03
rho_snow = 411.0
rho_water = 1000.0
c_water = 4.19E03
A_surf = 210.0
d_ins = 0.1
lam_i = 0.32

# Perform calculations
T_sol_air_vec = alpha * glob_solir_vec / h + air_temp_vec - T_cor_fact
Q_surf_vec = A_surf * lam_i / d_ins * T_sol_air_vec
f_srf_melt_vec = Q_surf_vec / (L_f * rho_snow)
hrly_srf_total_vec = f_srf_melt_vec * 3600.0

q_rain_vec = np.where(air_temp_vec > 0.0, prec_vec * rho_water * c_water * air_temp_vec / 3600.0, 0.0)
v_rain_vec = prec_vec * A_surf * rho_water * c_water * air_temp_vec / (L_f * rho_snow)
SMR_temp_vec = np.where(air_temp_vec > 0.0, hrly_srf_total_vec * rho_snow / A_surf, 0.0)
SMR_rain_vec = v_rain_vec * rho_snow / A_surf
SMR_total_vec = SMR_temp_vec + SMR_rain_vec
SMR_rainT_vec = np.cumsum(SMR_total_vec)
emp1_SMR_vec = -0.09 + 0.00014 * glob_solir_vec + 0.0575 * air_temp_vec + 0.0012 * air_temp_vec * air_vel_vec - 0.18 * air_temp_vec * d_ins
Psat_vec = Psat_WV(air_temp_vec + 273.15) / 10.0

# Extract the amount of RH precipitation column from the data
RH_perc_vec = rdata[:, 17].astype(float)
Pw_vec = Psat_vec * RH_perc_vec / 100.0  # kPa
w_vec = 2.16679 * Pw_vec * 1000 / (273.15 + air_temp_vec)  # kPa; 1000 to convert to kPa
emp2_SMR_vec = -0.97 - 0.097 * (d_ins * 100) + 0.164 * air_vel_vec + 0.00175 * glob_solir_vec + 0.102 * air_temp_vec + 0.192 * w_vec
emp1_SMR_wc_vec = np.where((emp1_SMR_vec < 0) | (air_temp_vec < 0), 0.0, emp1_SMR_vec)
emp2_SMR_wc_vec = np.where((emp2_SMR_vec < 0) | (air_temp_vec < 0), 0.0, emp2_SMR_vec)
emp1_SMR_wc_cs_vec = np.cumsum(emp1_SMR_wc_vec)
emp2_SMR_wc_cs_vec = np.cumsum(emp2_SMR_wc_vec)

# Display results using the custom print function
printNP(air_temp_vec, name="Air temperature (Celsius)")
printNP(q_rain_vec, name="Hourly melt rate from surface (m^3/h)")
printNP(SMR_rainT_vec, name="Rain and T cumulative (m^3/h)")

# Additional computations and vector operations

# Read another Parquet file using PyArrow
parquet_file = 'TsiTsoDATA.parquet'
table = pq.read_table(parquet_file)
data2 = table.to_pandas().to_numpy()
print(data2)

# Process additional columns
Tsi_vec_raw = data2[:, 0]
Tsi_vec = Tsi_vec_raw.astype(float)
printNP(Tsi_vec, name="Tsi (Celsius)")

Tso_vec_raw = data2[:, 1]
Tso_vec = Tso_vec_raw.astype(float)
printNP(Tso_vec, name="Tso (Celsius)")

# Calculate heat transfer coefficients
ho_vec = np.where(air_vel_vec <= 5.0, 6.0 + 4.0 * air_vel_vec, 7.41 * (air_vel_vec ** 0.78))
printNP(ho_vec, name="# Air velocity (with cond)")

# Internal surface heat transfer coefficient
h_i = 99.75  # W/m^2*K

# Heat flux calculations
qi_vec = (Tsi_vec - 0.0) * h_i  # W/m^2
printNP(qi_vec, name="Heat flux in (W/m^2)")

qo_vec = (T_sol_air_vec - Tso_vec) * ho_vec  # W/m^2
printNP(qo_vec, name="Heat flux out (W/m^2)")

# Phase change calculations
v_pc_vec = qi_vec / (L_f * rho_snow)  # m^3/(m^2*s)
printNP(v_pc_vec, name="Speed of phase change (m^3/(m^2*s))")

v_pc_hourly_vec = v_pc_vec * 3600  # m/h
printNP(v_pc_hourly_vec, name="Hourly speed of phase change (m^3/(m^2*h))")

hfmr_vec = np.where(air_temp_vec > 0, v_pc_hourly_vec * rho_snow, 0.0)
printNP(hfmr_vec, name="Hourly melt rate from solar heat flux")

hfmr_cumsum_vec = np.cumsum(hfmr_vec)
printNP(hfmr_cumsum_vec, name="Cumulative hourly melt rate from solar heat flux")

# Heat flux from rain and sun
rain_solar_hf_vec = q_rain_vec + qo_vec  # W/m^2
printNP(rain_solar_hf_vec, name="Heat flux from rain and sun (W/m^2)")

# Heat flux from wind, solar, and rain
wind_solar_rain_vec = rain_solar_hf_vec / (T_sol_air_vec - Tso_vec)
printNP(wind_solar_rain_vec, name="Heat flux from wind, solar and rain (W/m^2)")
