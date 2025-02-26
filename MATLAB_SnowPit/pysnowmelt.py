import time as t
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# Functions
def Geometry(Prop):
    a, b, alpha, beta, d, h = Prop[:6]

    # Precompute common terms
    tan_beta = np.tan(beta)
    tan_alpha = np.tan(alpha)
    sin_beta = np.sin(beta)
    sin_alpha = np.sin(alpha)

    # Volume calculations
    term_surf = (a - 2.0 * h / tan_beta) * (b - 2.0 * h / tan_beta)
    sqrt_surf = np.sqrt((a * b) * term_surf)
    Vol_surf = (h / 3.0) * ((a * b) + term_surf + sqrt_surf)

    term_below = (a - 2.0 * d / tan_alpha) * (b - 2.0 * d / tan_alpha)
    sqrt_below = np.sqrt((a * b) * term_below)
    Vol_below = (d / 3.0) * ((a * b) + term_below + sqrt_below)

    Vol_tot = Vol_below + Vol_surf

    # Area calculations
    a_surf = term_surf + (h / sin_beta) * (2.0 * a + 2.0 * b + 2.0 * \
            (a - 2.0 * h / tan_beta) + 2.0 * (b - 2.0 * h / tan_beta)) / 2.0
    a_below = term_below + (d / sin_alpha) * ((2.0 * a + 2.0 * b) + \
            2.0 * (a - 2.0 * d / tan_alpha) + 2.0 * (b - 2.0 * d / tan_alpha)) / 2.0

    return [Vol_surf, Vol_below, Vol_tot, a_surf, a_below]

#def birds_cycle(Table_weather):
#    Temp_in = np.zeros(365)
#    Temp_in[0] = 31  # Initial temperature
#    day = 0
#    for k in range(7):
#        end_day = day + 37
#
#        for i in range(day, end_day):
#            j = i + 1
#            Temp_in[j] = Temp_in[i] - 0.26
#
#        day = end_day + 1
#        rest = day + 14
#
#        for r in range(day, rest):
#            Temp_in[r] = Table_weather.iloc[r, 15]  # Max temp (r,9) or Heat Deg Days (r,15)???
#
#        Temp_in[rest] = 31
#
#        day = rest

def Melted_Volume(Temp, Rain, Prop, day):
    a, b, alpha, beta, d, h, Dw, Ds, Ls, cw, lambda_i, lambda_g, thickness = Prop

    # Precompute common terms
    tan_beta = np.tan(beta)
    tan_alpha = np.tan(alpha)
    sin_beta = np.sin(beta)
    sin_alpha = np.sin(alpha)

    # Surface area
    if h > 0:
        A = (a - 2.0 * h / tan_beta) * (b - 2.0 * h / tan_beta) + (h / sin_beta) * \
            (2.0 * a + 2.0 * b + 2.0 * (a - 2.0 * h / tan_beta) + 2.0 * (b - 2.0 * h / tan_beta)) / 2.0
    else:
        A = (a - 2.0 * abs(h) / tan_alpha) * (b - 2.0 * abs(h) / tan_alpha)

    # Ground area
    Ag = (a - 2.0 * d / tan_alpha) * (b - 2.0 * d / tan_alpha) + (d / sin_alpha) * \
        ((2.0 * a + 2.0 * b) + 2.0 * (a - 2.0 * d / tan_alpha) + 2.0 * (b - 2.0 * d / tan_alpha)) / 2.0

    # Surface melt
    Q_surf = A * (lambda_i / thickness) * Temp if Temp > 0 else 0
    V_surf = (Q_surf / (Ls * Ds)) * (3600.0 * 24.0)

    # Rain melt
    V_rain = ((Rain / 1000.0) * A * Dw * cw * Temp) / (Ls * Ds)

    # Ground melt
    Q_ground = lambda_g * Ag * 2.0
    V_ground = (Q_ground / (Ls * Ds)) * (3600.0 * 24.0)

    # Snowmelt for cooling
    cooling_days = {177: 0.01, 178: 0.015, 179: 0.020, 180: 0.025, 
                    181: 0.030, 182: 0.035, 183: 0.040, 184: 0.045}
    V_cool = cooling_days.get(day, 0) * 200

    # Total snowmelt
    V_melt = V_surf + V_rain + V_ground + V_cool
    V_melt0 = V_surf + V_rain + V_ground

    return V_melt, V_melt0, V_surf, V_rain, V_ground, V_cool

def findHeightNew(V, Prop):
    a = Prop[0]
    b = Prop[1]
    alpha = Prop[2]
    beta = Prop[3]
    d = Prop[4]

    Vol_below = (d / 3.0) * ((a * b) + (a - 2.0 * d / np.tan(alpha)) * (b - 2.0 * d / np.tan(alpha)) \
                + np.sqrt((a * b) * (a - 2.0 * d / np.tan(alpha)) * (b - 2.0 * d / np.tan(alpha))))

    if V > Vol_below:
        h_new = calculateHeight(V, Vol_below, a, b, beta)
    else:
        h_new = -calculateHeight(Vol_below - V, 0, a, b, alpha)

    return h_new

def calculateHeight(V, V_offset, a, b, angle):
    # Calculate the height of the snow pile using the volume and offset volume (V_offset) as input
    h_max = np.tan(angle) * a / 2.0 * 50.0  # 50 is a scaling factor 
    Vol_matrix = np.array([(i / 100.0 / 3.0) * ((a * b) + (a - 2.0 * i / 100.0 / np.tan(angle)) * \
        (b - 2.0 * i / 100.0 / np.tan(angle)) + np.sqrt((a * b) * (a - 2.0 * i / 100.0 / np.tan(angle)) *\
        (b - 2.0 * i / 100.0 / np.tan(angle)))) for i in range(1, int(h_max) + 1)])
    k = np.polyfit(Vol_matrix, np.arange(1, int(h_max) + 1), 1)
    h = (V - V_offset) / (k[0] * 100)
    return h

#############################################
################ MAIN CODE ##################
#############################################
def main():
    # Start timer
    start_time = t.time()

    print(
    " %%%%%%%%%%%%%%%%%%%%%%%%%%%% SNOW MELT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n",
     "% Code for simulating the snowmelt and volume variation of the SSS\n",
     "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    )
    
    # Initialize properties and load weather data
    print('\tInitializing properties and loading weather data...')
    a = 6.0  # Short side [m]
    b = 16.0  # Long side [m]
    alpha_degrees = 60.0  # Slope angle below surface
    beta_degrees = 45.0  # Slope angle above surface
    d = 2.0  # Pit depth [m]
    h = 2.0  # Pile height [m]
    beta = beta_degrees * np.pi / 180.0  # [rad]
    alpha = alpha_degrees * np.pi / 180.0  # [rad]
    Dw = 1000.0  # Density of water [kg m-3]
    Ds = 750.0  # Density of snow [kg m-3]
    Ls = 334000.0  # Latent heat of fusion [J kg-1]
    cw = 4185.0  # Heat capacity [J kg-1 K-1]
    lambda_i = 0.1  # Sawdust thermal conductivity [W m-1 K-1]
    lambda_g = 1.0  # Ground thermal conductivity [W m-1 K-1]
    thickness = 0.1  # Insulation layer thickness
    Prop = [a, b, alpha, beta, d, h, Dw, Ds, Ls, cw, lambda_i, lambda_g, thickness]
    Prop0 = Prop.copy()

    aa = Prop[0]
    bb = Prop[1]
    dd = Prop[4]
    hh = Prop[5]

    # Load weather data
    filename1 = 'Weather_Montreal_2018.xlsx'
    Table_weather = pd.read_excel(filename1)
    print('\tWeather data loaded successfully.')

    # Calculate initial geometry
    Geo = Geometry(Prop)
    Vi = Geo[0] + Geo[1]  # Initial volume of snow
    Vi0 = Geo[0] + Geo[1]  # Initial volume of snow
    mass = Vi * Ds

    Vf = np.zeros(6)
    Tair = Table_weather.iloc[:, 13].values  # Diurnal air mean temperature above freezing point [ºC]
    Tmax = Table_weather.iloc[:, 9].values
    P = Table_weather.iloc[:, 19].values  # Total Rain [mm]
    s = Tair.shape
    P[np.isnan(P)] = 0
    print('\tInitialization complete.')

    # Main loop with cooling
    print('\tStarting main loop with cooling...')
    num_of_days = s[0]
    Vf2 = np.zeros((len(Vf), num_of_days))  # Preallocating for speed
    Vnew2 = np.zeros(num_of_days)  # Preallocating for speed
    storeh = np.zeros(num_of_days)  # Preallocating for speed
    for i in range(s[0]):
        Temp = Tair[i]
        Rain = P[i]
        day = i + 1

        V_melt, V_melt0, V_surf, V_rain, V_ground, V_cool = Melted_Volume(Temp, Rain, Prop, day)

        Vf2[:, i] = [V_melt, V_melt0, V_surf, V_rain, V_ground, V_cool]
        Vnew2[i] = Vi - Vf2[0, i]  # NB! Vf2[0, i] = V_melt

        Vi = Vnew2[i]

        h_new = findHeightNew(Vnew2[i], Prop)  # More optimized function use
        storeh[i] = h_new

        Prop[5] = h_new
    print('\tMain loop with cooling completed.')

    # Main loop without cooling
    print('\tStarting main loop without cooling...')
    Vf0 = np.zeros((len(Vf), num_of_days))  # Preallocating for speed
    Vnew0 = np.zeros(num_of_days)  # Preallocating for speed
    storeh0 = np.zeros(num_of_days)  # Preallocating for speed
    for i in range(s[0]):
        Temp = Tair[i]
        Rain = P[i]
        day = i + 1

        V_melt, V_melt0, V_surf, V_rain, V_ground, V_cool = Melted_Volume(Temp, Rain, Prop0, day)

        Vf0[:, i] = [V_melt, V_melt0, V_surf, V_rain, V_ground, V_cool]
        Vnew0[i] = Vi0 - Vf0[1, i]  # NB! Vf0[1, i] = V_melt0

        Vi0 = Vnew0[i]

        h_new0 = findHeightNew(Vnew0[i], Prop)  # More optimized function use
        storeh0[i] = h_new0

        Prop0[5] = h_new0
    print('\tMain loop without cooling completed.')

    # End timer
    end_time = t.time()
    print(f"\nElapsed time for calculations: {end_time-start_time:.4e} sec.\n")

    # Plotting results
    print('Plotting results...')
    # Plot Figure 1
    f1, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    f1.suptitle('SNOWMELT VOLUME', fontsize=16)

    time = [datetime(2024, 1, 1) + timedelta(days=i) for i in range(num_of_days)]
    ax1.plot(time, storeh)
    ax1.set_ylabel('Height [m]')
    ax1.grid(True)
    ax1.legend(['HEIGHT'])

    data = Vf2.T
    ax2.plot(time, data[:, 0], color='#0072BD', label='TOTAL')
    ax2.plot(time, data[:, 2], color='#D95319', label='SURFACE')
    ax2.plot(time, data[:, 3], color='#EDB120', label='RAIN')
    ax2.plot(time, data[:, 4], color='#7E2F8E', label='GROUND')
    ax2.plot(time, data[:, 5], color='#77AC30', label='COOLING')
    ax2.set_ylabel('Volume [m^3]')
    ax2.grid(True)
    ax2.legend()

    # Plot Figure 2
    f2, ax3 = plt.subplots(figsize=(10, 6))
    f2.suptitle('PILE VOLUME', fontsize=16)

    i = np.where(Vnew2 < 0)[0]
    if len(i) == 0:
        ax3.plot(time, Vnew2, label='Volume WITH Cooling')
    else:
        Vnew2TOP = Vnew2.copy()
        Vnew0LOW = Vnew2.copy()
        Vnew2TOP[i[0]:] = np.nan
        Vnew0LOW[:i[0]] = np.nan
        ax3.plot(time, Vnew2TOP, label='Volume WITH Cooling')
        ax3.plot(time, Vnew0LOW, '--', color=ax3.get_lines()[-1].get_color())

    i = np.where(Vnew0 < 0)[0]
    if len(i) == 0:
        ax3.plot(time, Vnew0, label='Volume WITHOUT Cooling')
    else:
        Vnew0TOP = Vnew0.copy()
        Vnew0LOW = Vnew0.copy()
        Vnew0TOP[i[0]:] = np.nan
        Vnew0LOW[:i[0]] = np.nan
        ax3.plot(time, Vnew0TOP, '.', label='Volume WITHOUT Cooling')
        ax3.plot(time, Vnew0LOW, '--', color=ax3.get_lines()[-1].get_color())

    ax3.set_ylabel('Volume [m^3]')
    ax3.grid(True)
    ax3.legend()

    plt.show()
    print('Plotting completed.')

if __name__ == '__main__':
    main()
