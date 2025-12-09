import pandas as pd
import numpy as np
import numba as nb
import matplotlib.pyplot as plt

def Psat_WV(T_K):
    """
    Water vapour saturation pressure
    W. Wagner and A. Pruß:" The IAPWS Formulation 1995 for the
    Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use ",
    Journal of Physical and Chemical Reference Data, June 2002 ,Volume 31, Issue 2, pp.
    387535)

    Returns Saturation vapor pressure (hPa)
    ```
    # Example usage with a pandas DataFrame
    data = {'Temperature_K': [300, 310, 320, 330, 340]}  # Example temperatures in Kelvin
    df = pd.DataFrame(data)
    df['Saturation_Vapor_Pressure_hPa'] = Psat_WV_vectorized(df['Temperature_K'])
    print(df)
    ```
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
    
    x = Tc / T_K * (C1 * teta + C2 * np.power(teta, 1.5) + C3 * np.power(teta, 3) + C4 * np.power(teta, 3.5) + C5 * np.power(teta, 4) + C6 * np.power(teta, 7.5))
    
    x = np.exp(x) * Pc
    
    return x

@nb.jit(nopython=True)
def solve_tdma(a, b, c, d, n):
    """
    # Solve a tridiagonal system using the Thomas algorithm
    * a - lower diaganal
    * b - central
    * c - upper
    * d - coefficients
    * t - vector with unknowns
    * n - size of the system
    # Example usage
    ```
    aTest = np.array([0, 1, 1, 1])  # lower diagonal (a_1 to a_{n-1})
    bTest = np.array([4, 4, 4, 4])  # main diagonal (b_1 to b_n)
    cTest = np.array([1, 1, 1, 0])  # upper diagonal (c_1 to c_{n-1})
    dTest = np.array([5, 5, 5, 5])  # right-hand side
    nTest = len(dTest)
    solution = solve_tdma(aTest, bTest, cTest, dTest, nTest)
    print(solution)
    ```
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
    t_i = 0.0 # Inner temperature, C
    delta_db = d_ins # Layer thickness, m
    k_db = lam_i # Thermal conductivity
    a_db = D # Thermal diffusivity
    n_el = delta_db / dx # Number of numerical elements
    nr_hour = len(t_o) # Number of hours
    nodes = int(n_el+1)  # Number of nodes

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

            if k == 3500:
                k = k

            # Solve the tridiagonal system
            T_n = solve_tdma(a, b, c, d, nodes)

        # Store the result for this hour
        T_nh[:, h] = T_n

    return T_nh

# Constants
alpha = 0.8 # Solar light absorptivity
T_cor_fact = 4.0 # °C; correction factor for the temperature
h = 22.7 # W/(m^2K); heat transfer coefficient for external surface
L_f = 333.4E03  # J/kg; latent heat of fusion
rho_snow = 411.0  # kg/m^3; density of snow
rho_water = 1000  # kg/m^3; density of water
c_water = 4.19E03  # J/(kg*K); specific heat capacity of water
moist_cont = 50.0 # %; moisture content of snow
c_dry = 0.99E03 # J/(kg*K); specific heat capacity of dry snow
rho_dry = 100.0 # kg/m^3; density of dry snow
rho_wet = rho_dry + moist_cont/100.0 * 1000 # kg/m^3
lam = 1.0  # W/(mK); thermal conductivity of the ground
A_ground = 210.0  # m^2; area of the ground
d_ins = 0.1  # m; insulation layer thickness
lam_i = 0.32  # W/(mK); thermal conductivity for the insulating material
A_surf = 210.0  # m^2; surface area of the snow pile
h_i = 99.75  # W/m^2*K; heat transfer coefficient at the internal surface

# Load data
df = pd.read_parquet('SnowStorageDATA.parquet')
df2 = pd.read_parquet('SurfaceMeltRateDATA.parquet')
df3 = pd.read_parquet('TsiTsoDATA.parquet')

# Preprocess data
first_nan_index = df['YEAR.1'].isna().idxmax()
df_wo_nan = df.loc[:first_nan_index - 1]

period_start_in = df[df['time.1'] == '2023-04-01T00:00'].index[0]
period_end_in = df[df['time.1'] == '2023-08-31T23:00'].index[0]
df_wo_nan_period = df_wo_nan.loc[period_start_in:period_end_in].reset_index(drop=True)

# Extract relevant columns
df_air_vel = df_wo_nan_period.iloc[:, 15]
df_air_temp = df_wo_nan_period.iloc[:, 7]
df_am_of_perc = df_wo_nan_period.iloc[:, 5]  # m/h
df_Gsi = df_wo_nan_period['ALLSKY_SFC_SW_DWN']
df_RH_perc = df_wo_nan_period.iloc[:, 17]
df_Tsi_from_file = df3['Tsi_C']
df_Tso_from_file = df3['Tso_C']

# Calculations
Q_ground = lam * A_ground * (6.0 - 0.0) / 1.0  # W
Q_melt = Q_ground / (L_f * rho_snow)  # m^3/s
Q_melt_total = Q_melt * 3600  # m^3/h

df_T_sol_air = alpha * df_Gsi / h + df_air_temp - T_cor_fact
df_Q = A_surf * lam_i / d_ins * (df_T_sol_air - 0.0)  # W
c_wet = (1.0 - moist_cont/100.0)*c_dry + moist_cont/100.0*c_water
print(f"Specific heat capacity (wet): {c_wet/1000.0:.4} kJ/(kg*K).")

D = lam_i/(c_wet * rho_wet) # m^2/s
print(f"Thermal diffusivity of the insulating material is {D:.4e} m^2/s.")

df_f_melt_surf = df_Q / (L_f * rho_snow)  # m^3/s
df_f_melt_surf_total = df_f_melt_surf * 3600  # m^3/h

df_q_rain = pd.Series(0.0, index=df_air_temp.index)
positive_temp_mask = df_air_temp > 0
df_q_rain[positive_temp_mask] = df_am_of_perc[positive_temp_mask] * rho_water * c_water * df_air_temp[positive_temp_mask] / 3600.0

df_v_rain = df_am_of_perc * A_surf * rho_water * c_water * df_air_temp / (L_f * rho_snow)  # m^3/h
df_SMR_temp = pd.Series(0.0, index=df_air_temp.index)
smr_temp_mask = df_air_temp > 0.0
df_SMR_temp[smr_temp_mask] = df_f_melt_surf_total[smr_temp_mask] * rho_snow / A_surf  # kg/(m^2*h)

df_SMR_rain = df_v_rain * rho_snow / A_surf  # kg/(m^2*h)
df_SMR_total = df_SMR_temp + df_SMR_rain
df_SMR_total_rain_cs = df_SMR_total.cumsum()

# Calculate the saturation vapor pressure at the air temperature
df_Psat = Psat_WV(df_air_temp + 273.15)/10.0 # hPa; 100/1000 to convert to hPa
df_Pw = df_Psat*df_RH_perc/100.0 # kPa; # Water steam pressure

# Absolute humidity
df_w = 2.16679*df_Pw*1000/(273.15+df_air_temp) # g/m^3; 1000 to convert to kPa

# Air velocity with conditions
df_ho = pd.Series(np.where(
        df_air_vel <= 5.0,
        6.0 + 4.0*df_air_vel,
        7.41 * (df_air_vel**0.78)
    ))

t_o = df_T_sol_air.to_numpy()
h_o = df_ho.to_numpy()
t_o_range = transient1D(t_o, h_o, d_ins, lam_i, D)

df_Tsi = pd.Series(t_o_range[-1, :])
df_Tso = pd.Series(t_o_range[ 0, :])

df_qi = (df_Tsi - 0.0)* h_i # W/m^2; # Heat flux in
df_qo = (df_T_sol_air - df_Tso) * df_ho # W/m^2; # Heat flux out
df_v_pc = df_qi/(L_f * rho_snow) # m^3/(m^2*s)
df_v_pc_hourly = df_v_pc * 3600
df_hfmr = pd.Series(np.where(df_air_temp > 0,  df_v_pc_hourly * rho_snow, 0.0))
df_hfmr_cs = df_hfmr.cumsum()
df_rain_solar_hf = df_q_rain + df_qo    # W/m^2
df_wind_solar_rain = df_rain_solar_hf/ (df_T_sol_air - df_Tso)

# Empirical models
df_SMR_emp1 = -0.09 + 0.00014 * df_Gsi + 0.0575*df_air_temp + 0.0012*df_air_temp*df_air_vel - 0.18*df_air_temp*d_ins # kg/m2/h
df_SMR_emp1_cond = pd.Series(np.where((df_SMR_emp1 < 0) | (df_air_temp < 0), 0, df_SMR_emp1))
df_SMR_emp1_cond_cs = df_SMR_emp1_cond.cumsum()

df_SMR_emp2 = -0.97 - 0.097*(d_ins*100) + 0.164*df_air_vel + 0.00175*df_Gsi + 0.102*df_air_temp + 0.192*df_w # kg/m2/h
df_SMR_emp2_cond = pd.Series(np.where((df_SMR_emp2 < 0) | (df_air_temp < 0), 0, df_SMR_emp2))
df_SMR_emp2_cond_cs = df_SMR_emp2_cond.cumsum()

# Printing
print("Air velocity:\n",df_air_vel)
print("\nAir temperature:\n", df_air_temp)
print("\nAmount of percipitation:\n", df_am_of_perc) 
print("\nGlobalr solar irradiance W/m2:\n", df_Gsi)
print("\nRelative humidity of percipitation:\n", df_RH_perc)
print("\nAir velocity with conditions:\n", df_ho)
print("\nInternal tempterature of snow:\n", df_Tsi)
print("\nOuter tempterature of snow:\n", df_Tso)
print("\nInternal layer heat flux  (W/m^2):\n", df_qi)
print("\nOuter layer heat flux put (W/m^2):\n", df_qo)
#print(df_v_pc)
print("\nHourly speed of phase change (m^3/(m^2*h)):\n", df_v_pc_hourly)
#print(df_hfmr)
print("\nCumulative hourly melt rate from solar heat flux:\n", df_hfmr_cs)
print("\nHeat flux from rain and sun (W/m^2):\n", df_rain_solar_hf)
print("\nHeat flux from wind, solar and rain (W/m^2)\n:",df_wind_solar_rain)
#print(df_SMR_emp1)
#print(df_SMR_emp1_cond)
print("\nCumulative Empirical 1:\n", df_SMR_emp1_cond_cs)
#print(df_SMR_emp2)
#print(df_SMR_emp2_cond)
print("\nCumulative Empirical 2:\n", df_SMR_emp2_cond_cs)


## Plotting
plt.figure(figsize=(8, 6))
plt.scatter(df_air_temp, df_SMR_emp2_cond, label='Surface 2')
plt.scatter(df_air_temp, df_SMR_total, label='Total Rain')
plt.scatter(df_air_temp, df_hfmr, label='Total Solar')
plt.scatter(df_air_temp, df_SMR_emp1_cond, label='Surface 1')
plt.title('Scatter Plot from Two DataFrames')
plt.xlabel('Air temperature C')
plt.ylabel('Melt Rate')
plt.legend()
plt.grid()
plt.show()

# Time series plot
df_time = pd.to_datetime(df_wo_nan_period['time']).astype('int64') / 1e9
fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(df_time, df_T_sol_air, label='Solar - Air (Left Axis)')
ax1.plot(df_time, df_Tso, color='red', label='Outer Surface Temp (Left Axis)')
ax1.set_xlabel('Time')
ax1.set_ylabel('Temperature C')
ax1.set_title('Change of temperature over time')
ax1.legend(loc='upper left')

ax2 = ax1.twinx()
ax2.plot(df_time, df_Tsi, label='Inner Surface Temp (Right Axis)', color='orange')
#ax2.set_ylim(-0.1, 0.8)
ax2.set_ylim(-0.5, 2.0)
ax2.set_ylabel('Solar In (Temperature)')
ax2.legend(loc='lower right')
plt.show()

# Cumulative melt rate plot
plt.figure(figsize=(8, 6))
plt.plot(df_time, df_SMR_emp1_cond_cs, label='Empirical 1')
plt.plot(df_time, df_SMR_emp2_cond_cs, label='Empirical 2')
plt.plot(df_time, df_SMR_total_rain_cs, label='Total Rain')
plt.plot(df_time, df_hfmr_cs, label='Total Solar')
plt.title('Cumulative melt rate over time')
plt.xlabel('Time')
plt.ylabel('Surface Melt Rate, kg/m^2')
plt.legend()
plt.grid()
plt.show()