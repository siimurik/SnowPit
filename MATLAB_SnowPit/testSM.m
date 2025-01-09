clc
clear

a = 6; % Short side [m]
b = 16; % Long side [m]
alpha_degrees = 60; % Slope angle below surface
beta_degrees = 45; % Slope angle above surface
d = 2; % Pit depth
h = 2; % Pile height
beta = beta_degrees * pi / 180;
alpha = alpha_degrees * pi / 180;
Dw = 1000; % Density of water [kg m-3]
Ds = 750; % Density of snow [kg m-3]
Ls = 334000; % Latent heat of fusion [J kg-1]
cw = 4185; % Heat capacity [J kg-1 K-1]
lambda_i = 0.1; % Sawdust thermal conductivity [W m-1 K-1]
lambda_g = 1; % Ground thermal conductivity [W m-1 K-1]
thickness = 0.1; % Insulation layer thickness
Prop = [a b alpha beta d h Dw Ds Ls cw lambda_i lambda_g thickness];
Prop0 = Prop;

aa = Prop(1);
bb = Prop(2);
dd = Prop(5);
hh = Prop(6);

%% https://climate.weather.gc.ca/historical_data/search_historic_data_e.html
filename1 = fullfile('Weather_Montreal_2018.xlsx');
Table_weather = readmatrix(filename1); % Open Excel file

Geo = Geometry(Prop);
Vi  = [Geo(1) + Geo(2)]; % Initial volume of snow
Vi0 = [Geo(1) + Geo(2)]; % Initial volume of snow
mass = Vi * Ds;

Vf = [0 0 0 0 0 0]';
Tair = Table_weather(:, 14); % Diurnal air mean temperature above freezing point [ºC]
Tmax = Table_weather(:, 10);
P = Table_weather(:, 20); % Total Rain [mm]
s = size(Tair);
P(isnan(P)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
birds_cycle(Table_weather);

% LOOP VARYING H WITH COOLING
for i = 1:s(1)
    Temp = Tair(i);
    Rain = P(i);
    day = i;

    [V_melt, V_melt0, V_surf, V_rain, V_ground, V_cool] = Melted_Volume(Temp, Rain, Prop, day);

    Vf2(:, i) = [V_melt V_melt0 V_surf V_rain V_ground V_cool]'; % New Volume of the pile
    Vnew2(i) = Vi - Vf2(1, i);

    Vi = Vnew2(1, i);

    h_new = findHeight(Vnew2(i), Prop); % New height of the pile
    storeh(i) = h_new;

    Prop(6) = h_new;
end


%%
% FUNCTION TO CALCULATE DIMENSIONS OF THE PILE
function Geo = Geometry(Prop)
    a = Prop(1);
    b = Prop(2);
    alpha = Prop(3);
    beta = Prop(4);
    d = Prop(5);
    h = Prop(6);

    % VOLUME
    Vol_surf = (h/3) * ( (a*b) + (a - 2*h/tan(beta)) * (b - 2*h/tan(beta)) + sqrt((a*b)*(a - 2*h/tan(beta)) * (b - 2*h/tan(beta))));
    Vol_below = (d/3) * ( (a*b) + (a - 2*d/tan(alpha)) * (b -2*d/tan(alpha)) + sqrt((a*b)*(a - 2*d/tan(alpha)) * (b -2*d/tan(alpha))));
    Vol_tot = Vol_below + Vol_surf;

    % AREA
    a_surf = (a - 2*h/tan(beta)) * (b - 2*h/tan(beta)) + (h/sin(beta))*(2*a + 2*b + 2*(a - 2*h/tan(beta)) + 2*(b -2*h/tan(beta)))/2;
    a_below = (a - 2*d/tan(alpha)) * (b - 2*d/tan(alpha)) +(d/sin(alpha))*((2*a + 2*b) + 2*(a - 2*d/tan(alpha)) + 2*(b -2*d/tan(alpha)))/2;

    Geo = [Vol_surf Vol_below Vol_tot a_surf a_below];
end

function birds_cycle(Table_weather)
    Temp_in = [31];
    day = 1;
    for k = 1:7
    end_day = day + 37;

    for i = day:end_day
        j = i + 1;
        Temp_in(j) = Temp_in(i) - 0.26;
    end

    day = end_day + 1;
    rest = day + 14;

    for r = day:rest
        Temp_in(r) = Table_weather(r,16);   % Max temp (r,10) or Heat Deg Days (r,16)???
    end

    Temp_in(rest) = 31;

    day = rest;
    end
end

% FUNCTION TO CALCULATE MELTED VOLUME
function [V_melt, V_melt0, V_surf, V_rain, V_ground, V_cool] = Melted_Volume(Temp,Rain,Prop,day)
    a = Prop(1);
    b = Prop(2);
    alpha = Prop(3);
    beta = Prop(4);
    d = Prop(5);
    h = Prop(6);
    Dw = Prop(7);
    Ds = Prop(8);
    Ls = Prop(9);
    cw = Prop(10);
    lambda_i = Prop(11);
    lambda_g = Prop(12);
    thickness = Prop(13);
   
    % SURFACE MELT
    % Insulated with wood chips
    if h > 0
        A = (a - 2*h/tan(beta)) * (b - 2*h/tan(beta)) + (h/sin(beta))*((2*a   + 2*b) + 2*(a - 2*h/tan(beta)) + 2*(b - 2*h/tan(beta)))/2;
    else
        A = (a - 2*abs(h)/tan(alpha)) * (b - 2*abs(h)/tan(alpha));
    end
   
    % Ag = (d/3) * ( (a*b) + (a - 2*d/tan(alpha)) * (b - 2*d/tan(alpha))+ sqrt((a*b)*(a - 2*d/tan(alpha)) * (b - 2*d/tan(alpha))));
    Ag = (a - 2*d/tan(alpha)) * (b - 2*d/tan(alpha)) +   (d/sin(alpha))*((2*a + 2*b) + 2*(a - 2*d/tan(alpha)) + 2*(b -   2*d/tan(alpha)))/2;
   
    if Temp <= 0
        Q_surf = 0;
    else
        Q_surf = A * (lambda_i / thickness) * (Temp); % Thermal flow [W]
    end
   
    v_surf = Q_surf / (Ls * Ds); % [m3 s-1]
    V_surf = v_surf * (3600 * 24); % [m3 day-1]
   
    % RAIN MELT
    V_rain = ((Rain/1000) * A * Dw * cw * (Temp)) / (Ls * Ds); % P[m],T[K], cs [J kg-1 K-1]
   
    % GROUND MELT
    Q_ground = lambda_g * Ag * (2); % Thermal flow [W]
    v_ground = Q_ground / (Ls * Ds); % [m3 s-1]
    V_ground = v_ground * (3600 * 24); % [m3 day-1]
   
    % SNOWMELT FOR COOLING
   
    if day == 177 || day == 225
        V_cool = 0.01 * 200;
    elseif day == 178 || day == 226
        V_cool = 0.015 * 200;
    elseif day == 179 || day == 227
        V_cool = 0.020 * 200;
    elseif day == 180 || day == 228
        V_cool = 0.025 * 200;
    elseif day == 181 || day == 229
        V_cool = 0.030 * 200;
    elseif day == 182 || day == 230
        V_cool = 0.035 * 200;
    elseif day == 183 || day == 231
        V_cool = 0.040 * 200;
    elseif day == 184 || day == 232
        V_cool = 0.045 * 200;
    %elseif day == 185 | day == 233
        %V_cool = 0.050 * 200;
    else
        V_cool = 0;
    end


    % TOTAL SNOWMELT

    V_melt = V_surf + V_rain + V_ground + V_cool;
    V_melt0 = V_surf + V_rain + V_ground;

end

% FUNCTION TO CALCULATE THE HEIGHT OF THE PILE FOR A GIVEN VOLUME
function h_new = findHeight(V,Prop)
    a = Prop(1);
    b = Prop(2);
    alpha = Prop(3);
    beta = Prop(4);
    d = Prop(5);

    Vol_below = (d/3) * ( (a*b) + (a - 2*d/tan(alpha)) * (b - 2*d/tan(alpha)) + sqrt((a*b)*(a - 2*d/tan(alpha)) * (b - 2*d/tan(alpha))));

    if V > Vol_below

        Vol_matrix = [0];
        h_max = tan(beta) * a / 2;
        h_max = h_max * 50;
        for i = 1:h_max
            h = i/100;
            Vol = (h/3) * ( (a*b) + (a - 2*h/tan(beta)) * (b -2*h/tan(beta)) + sqrt((a*b)*(a - 2*h/tan(beta)) * (b - 2*h/tan(beta))));
            Vol_matrix(i) = Vol;
        end
        [r,m,b] = fitlm(1:h_max,Vol_matrix); 
        h_new = (V - Vol_below) / (m*100);

    else

        Vol_matrix = [0];
        d_max = tan(alpha) * a / 2;
        d_max = d_max * 50;
        for i = 1:d_max
            d = i/100;
            Vol = (d/3) * ( (a*b) + (a - 2*d/tan(alpha)) * (b -2*d/tan(alpha)) + sqrt((a*b)*(a - 2*d/tan(alpha)) * (b -2*d/tan(alpha))));
            Vol = Vol_below - Vol;
        Vol_matrix(i) = Vol;
        end

        [r,m,b] = fitlm(1:d_max, Vol_matrix);
        h_new = -(V-Vol_below)/(m*100);

    end
end
