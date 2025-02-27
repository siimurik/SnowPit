# pip install numpy pyarrow numba
import numpy as np
import pyarrow.parquet as pq
import numba as nb

def printNP(arr, name=None):
    """
    Prints a NumPy array in the usual format with optional name, length, and dtype.
    """
    array_str = np.array2string(arr, separator=' ', threshold=6)
    length = arr.shape[0]
    dtype = arr.dtype
    if name:
        print(f"\n{array_str}\nName: {name}, Length: {length}, dtype: {dtype}")
    else:
        print(f"\n{array_str}\nLength: {length}, dtype: {dtype}")

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
    
    x = Tc / T_K * (C1*teta + C2*teta**1.5 + C3*teta**3 + C4*teta**3.5 + C5*teta**4 + C6*teta**7.5)
    
    x = np.exp(x) * Pc
    
    return x

@nb.jit(nopython=True)
def solve_tdma(a, b, c, d, n):
    """
    Solves a tridiagonal system of equations using the Thomas algorithm.

    Parameters:
        a (np.array): Lower diagonal of the matrix.
        b (np.array): Main diagonal of the matrix.
        c (np.array): Upper diagonal of the matrix.
        d (np.array): Right-hand side vector.
        n (int): Size of the system.

    Returns:
        np.array: Solution vector.
    """
    n = int(n)
    c_prime = np.zeros(n)
    d_prime = np.zeros(n)
    x = np.zeros(n)

    # Forward sweep
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]

    for i in range(1, n):
        c_prime[i] = c[i] / (b[i] - a[i] * c_prime[i - 1])
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / (b[i] - a[i] * c_prime[i - 1])

    # Back substitution
    x[-1] = d_prime[-1]

    for i in range(n-2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]

    return x

@nb.jit(nopython=True)
def transient1D(t_o, h_o, d_ins, lam_i, D, dx=0.005, dt=10.0, h_i=99.75):
    """
    Temperature distribution, transient 1D, BC 3d, implicit method
    """
    t_i = 0.0 # Inner temperature, C
    delta_db = d_ins # Layer thickness, m
    k_db = lam_i # Thermal conductivity
    a_db = D # Thermal diffusivity
    n_el = delta_db / dx # Number of numerical elements
    nr_hour = len(t_o) # Number of hours
    nodes = int(n_el + 1)  # Number of nodes

    # Initialize temperature distribution
    T_n = np.zeros(nodes)
    T_nh = np.zeros((nodes, nr_hour))

    # Initialize arrays for Tri-Diagonal Matrix Algorithm (TDMA)
    a = np.zeros(nodes)
    b = np.zeros(nodes)
    c = np.zeros(nodes)
    d = np.zeros(nodes)

    nh = int(3600 / dt) # Number of time intervals in one hour

    # Main loop over hours
    for h in range(nr_hour):
        dFo = a_db * dt / dx**2
        dBio_i = h_i * dx / k_db
        dBio_o = h_o[h] * dx / k_db

        # Time steps within one hour
        for k in range(nh):
            # Set up the tridiagonal system
            b[0] = (1.0 + 2.0*dFo + 2.0*dFo*dBio_o)
            c[0] = -2.0 * dFo
            d[0] = T_n[0] + 2.0*dFo*dBio_o*t_o[h]

            for j in range(1, nodes-1):
                a[j] = -dFo
                b[j] = (1.0 + 2.0*dFo)
                c[j] = -dFo
                d[j] = T_n[j]

            a[-1] = -2.0 * dFo
            b[-1] = (1.0 + 2.0*dFo + 2.0*dFo*dBio_i)
            d[-1] = T_n[-1] + 2.0*dFo*dBio_i*t_i

            # if k == 3500:   # ???
            #     k = k       # ???

            # Solve the tridiagonal system
            T_n = solve_tdma(a, b, c, d, nodes)

        # Store the result for this hour
        T_nh[:, h] = T_n

    return T_nh


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
emp1_SMR_vec = -0.09 + 0.00014*glob_solir_vec + 0.0575*air_temp_vec + 0.0012*air_temp_vec*air_vel_vec - 0.18*air_temp_vec*d_ins
Psat_vec = Psat_WV(air_temp_vec + 273.15) / 10.0

# Extract the amount of RH precipitation column from the data
RH_perc_vec = rdata[:, 17].astype(float)
Pw_vec = Psat_vec * RH_perc_vec / 100.0  # kPa
w_vec = 2.16679 * Pw_vec * 1000 / (273.15 + air_temp_vec)  # kPa; 1000 to convert to kPa
emp2_SMR_vec = -0.97 - 0.097*(d_ins*100) + 0.164*air_vel_vec + 0.00175*glob_solir_vec + 0.102*air_temp_vec + 0.192*w_vec
emp1_SMR_wc_vec = np.where((emp1_SMR_vec < 0) | (air_temp_vec < 0), 0.0, emp1_SMR_vec)
emp2_SMR_wc_vec = np.where((emp2_SMR_vec < 0) | (air_temp_vec < 0), 0.0, emp2_SMR_vec)
emp1_SMR_wc_cs_vec = np.cumsum(emp1_SMR_wc_vec)
emp2_SMR_wc_cs_vec = np.cumsum(emp2_SMR_wc_vec)

# Display results using the custom print function
printNP(air_temp_vec, name="Air temperature (Celsius)")
printNP(q_rain_vec, name="Hourly melt rate from surface (m^3/h)")
printNP(SMR_rainT_vec, name="Rain and T cumulative (m^3/h)")

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
ho_vec = np.where(air_vel_vec <= 5.0, 6.0 + 4.0 * air_vel_vec, 7.41 * (air_vel_vec**0.78))
printNP(ho_vec, name="# Air velocity (with cond)")

# More constants for transient 1D solver
h_i = 99.75  # W/m^2*K
c_wet = 2.59E03 # J/(kg*K)
rho_dry = 100.0 # kg/m^3
moist_cont = 50.0 # %
rho_wet = rho_dry + moist_cont/100.0*1000 # kg/m^3
c_dry = 0.99E03 # J/(kg*K)
c_water = 4.19E03 # J/(kg * K)

c_wet = (1.0 - moist_cont/100.0)*c_dry + moist_cont/100.0*c_water
print(f"Specific heat capacity (wet): {c_wet/1000.0:.4} kJ/(kg*K).")

D = lam_i/(c_wet * rho_wet) # m^2/s
print(f"Thermal diffusivity of the insulating material is {D:.4e} m^2/s.")

# Snow outer and inner layer temperatures
t_o = T_sol_air_vec
h_o = ho_vec

t_o_range = transient1D(t_o, h_o, d_ins, lam_i, D, dx=0.005, dt=10.0, h_i=h_i)
printNP(t_o_range[0, :], name="Outer temperature of snow (°C)")  
printNP(t_o_range[-1, :], name="Internal temperature of snow (°C)")  

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
