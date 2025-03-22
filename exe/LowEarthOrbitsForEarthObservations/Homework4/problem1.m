%% Preamble
clc; clear all; close all;

%% Provided constants
earth_semimajor_axis = 6378136.3;
earth_gravitational_parameter = 3.986004415e14;
earth_angular_velocity = 7.292115e-5;
earth_oblateness_parameter = 0.001082;

%% Orbit parameters
orbit_semimajor_axes = [
    earth_semimajor_axis + 717e3;
    earth_semimajor_axis + 1336e3;
    earth_semimajor_axis + 490e3;
    earth_semimajor_axis + 782e3;
    earth_semimajor_axis + 475e3;
    earth_semimajor_axis + 812e3;
];

orbit_eccentricities = [
    0.0003565;
    0.0007665;
    0.00005855;
    0.0031394;
    0.0004512;
    0.0205681;
];

orbit_inclinations = deg2rad([
    92.0194;
    66.0425;
    88.9691;
    98.7932;
    97.4358;
    49.8242;
]);

%% Check which satellites are in sun-synchronous orbits
% Rename
Re = earth_semimajor_axis;
mu = earth_gravitational_parameter;
J2 = earth_oblateness_parameter;
a = orbit_semimajor_axes;
e = orbit_eccentricities;
I = orbit_inclinations;

sun_synchronous_precession_in_RAAN = deg2rad(0.00416667);

n_bar = sqrt(mu ./ a.^3);

Omega_bar_dot = ...
    -(3/2) .* n_bar .* (Re ./ a).^2 .* J2 .* cos(I) ./ (1 - e.^2).^2;

omega_bar_dot = ...
    -(3/4) .* n_bar .* (Re ./ a).^2 .* J2 .* (1 - 5 .* cos(I).^2) ./ (1 - e.^2).^(2);

M_bar_dot = ...
    n_bar .* (1 - (3/4) .* (Re ./ a).^2 .* J2 .* (1 - 3 .* cos(I).^2) ./ (1 - e.^2).^(3/2));

Omega_s_dot     = 1.99096871e-7;
theta_g_dot = deg2rad(360 / 86164.0905);

keplerian_period      = 2 * pi * sqrt(a.^3 ./ mu);
anomalistic_period    = 2 * pi ./ M_bar_dot;
draconitic_period     = 2 * pi ./ (M_bar_dot + omega_bar_dot);
nodal_day             = 2 * pi ./ (Omega_bar_dot - theta_g_dot);
sun_cycle_period      = 2 * pi ./ (Omega_bar_dot - Omega_s_dot);
sun_cycle_period_days = sun_cycle_period ./ (60 * 60 * 24);

%% Investigate the Beta prime angle
max_non_sun_synchronous_sun_cycle_period_days = max(abs(sun_cycle_period_days(abs(sun_cycle_period_days) < 1000)));

max_time_days = 4 * max_non_sun_synchronous_sun_cycle_period_days;

time_days = 0:1:round(max_time_days);
time_seconds = time_days * (60 * 60 * 24);

% Satellite longitude of ascending node
Omega_0 = 0;
Omega_hist = Omega_0 + Omega_bar_dot .* reshape(time_seconds, 1, []);

% Sun right ascension
Omega_s_0 = 0;
Omega_s_hist = Omega_s_0 + Omega_s_dot .* reshape(time_seconds, 1, []);

% Ecliptic tilt
epsilon_0 = deg2rad(23.439);

num_sats = size(Omega_hist, 1);
num_time = numel(time_seconds);

beta_prime = zeros(num_sats, num_time);

for sat_idx = 1:num_sats
    for time_idx = 1:num_time

        % Compute e_hat_h
        Rx = ValladoROT1(rad2deg(orbit_inclinations(sat_idx)));
        Rz = ValladoROT3(rad2deg(Omega_hist(sat_idx, time_idx)));

        e_hat_h = transpose(Rz) * transpose(Rx) * [0; 0; 1];

        % Compute e_hat_sun
        Rx_sun = ValladoROT1(rad2deg(epsilon_0));
        Rz_sun = ValladoROT3(rad2deg(Omega_s_hist(time_idx)));

        e_hat_sun = transpose(Rx_sun) * transpose(Rz_sun) * [1; 0; 0];

        % Compute beta prime angle
        beta_prime(sat_idx, time_idx) = pi/2 - acos(dot(e_hat_h, e_hat_sun));
    end
end


%% Analysis
[~, locs_1] = findpeaks(beta_prime(1, :));
[~, locs_2] = findpeaks(beta_prime(2, :));
[~, locs_3] = findpeaks(beta_prime(3, :));
[~, locs_6] = findpeaks(beta_prime(6, :));

time_of_peaks_1 = time_days(locs_1);
time_of_peaks_2 = time_days(locs_2);
time_of_peaks_3 = time_days(locs_3);
time_of_peaks_6 = time_days(locs_6);

period_1 = diff(time_of_peaks_1);
period_2 = diff(time_of_peaks_2);
period_3 = diff(time_of_peaks_3);
period_6 = diff(time_of_peaks_6);



%% Plotting
figure('Position', [100, 100, 1500, 750]);
tlt = tiledlayout(4, 1);
title(tlt, "\beta ' Angles for Cryosat-2, TOPEX/Poseidon, GRACE-FO-1, and Starlette", ...
    "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");

nexttile(1);
h1 = plot(time_days, rad2deg(beta_prime(1, :)), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Time [days]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("\beta '", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylim([-100, 100]);

nexttile(2);
h1 = plot(time_days, rad2deg(beta_prime(2, :)), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Time [days]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("\beta '", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylim([-100, 100]);

nexttile(3);
h1 = plot(time_days, rad2deg(beta_prime(3, :)), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Time [days]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("\beta '", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylim([-100, 100]);

nexttile(4);
h1 = plot(time_days, rad2deg(beta_prime(6, :)), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Time [days]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("\beta '", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylim([-100, 100]);



