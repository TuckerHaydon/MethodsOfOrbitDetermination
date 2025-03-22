%% Preamble
clc; clear all; close all;

%% Constants
J2 = 0.001082;
ae = 6378136.3;
mu = 3.986004415e14;
sidereal_day_sec = 86164.1;
sidereal_rate_deg = 360 / sidereal_day_sec;

%% Given
TRMM_eccentricity        = 0;
TRMM_inclination         = deg2rad(35);
TRMM_semimajor_axis      = ae + 350e3;
TRMM_time_ascending_node = datetime(1999, 01, 21, 20, 43, 47);
TRMM_longitude           = deg2rad(-5.157);
TRMM_longitude_deg       = rad2deg(TRMM_longitude);

% Sun synchronous orbit
% Resurs_eccentricity     = 

%% Calculations
TRMM_mean_motion  = sqrt(mu / TRMM_semimajor_axis^3);

% Rename
n_bar = TRMM_mean_motion;
Re = ae;
a = TRMM_semimajor_axis;
I = TRMM_inclination;
e = TRMM_eccentricity;

omega_bar_dot = ...
    -(3/4) .* n_bar .* (Re ./ a).^2 .* J2 .* (1 - 5 .* cos(I).^2) ./ (1 - e.^2).^(2);

M_bar_dot = ...
    n_bar .* (1 - (3/4) .* (Re ./ a).^2 .* J2 .* (1 - 3 .* cos(I).^2) ./ (1 - e.^2).^(3/2));


draconitic_period = 2 * pi ./ (M_bar_dot + omega_bar_dot);

Omega_bar_dot = -(3/2) .* n_bar .* (Re ./ a).^2 .* J2 .* cos(I) ./ (1 - e.^2).^2;
Omega_bar_dot_deg = rad2deg(Omega_bar_dot);

local_longitude_regression_deg = Omega_bar_dot_deg - sidereal_rate_deg;

N = 0:1:6000;

future_crossing_times = TRMM_time_ascending_node + seconds(N * draconitic_period);

time_mask = year(future_crossing_times) == 1999;

N = N(time_mask);
future_crossing_times = future_crossing_times(time_mask);

future_crossing_longitudes_deg = TRMM_longitude_deg + wrapTo180(local_longitude_regression_deg * N * draconitic_period);
future_crossing_times_offset_hours = future_crossing_longitudes_deg / 15;

future_crossing_times_LMT = future_crossing_times + hours(future_crossing_times_offset_hours);

LMT_hours = hour(future_crossing_times_LMT) + minute(future_crossing_times_LMT) / 60;

LMT_mask = abs(LMT_hours - 22.33333333) < (1/60);

figure('Position', [100, 100, 1000, 750]); 
hold on;
plot(future_crossing_times, LMT_hours, "k-", "LineWidth", 3.0);
xlims = xlim();
plot(xlims, [22.333, 22.333], "r--", "LineWidth", 3.0);
grid on;
grid minor;
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Node Crossing Number [Unitless]", "FontSize", 14, "FontWeight", "bold", "FontName", "Times New Roman");
ylabel("Local Mean Time [hours]", "FontSize", 14, "FontWeight", "bold", "FontName", "Times New Roman");
% ylim([-10, 30]);
title("TRMM Ascending Node Crossing Times in Local Mean Time", "FontSize", 18, "FontWeight", "bold", "FontName", "Times New Roman")



