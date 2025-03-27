%% Preamble
clc; clear all; close all;

%% Provided values
% Time
initial_utc_year    = 2018;
initial_utc_month   = 02;
initial_utc_day     = 1;
initial_utc_hour    = 05;
initial_utc_minute  = 00;
initial_utc_seconds = 0;
initial_utc_time    = datetime(...
    initial_utc_year, ...
    initial_utc_month, ...
    initial_utc_day, ...
    initial_utc_hour, ...
    initial_utc_minute, ...
    initial_utc_seconds);

eop_file = "finals.all.iau1980.txt";
[...      
        ut1_utc_sec, ...
        polar_motion_deg, ...
        length_of_day_earth_orientation_parameter_sec, ...
        longitude_earth_orientation_parameter_deg, ...
        obliquity_earth_orientation_parameter_deg] = ...
    ParseEarthOrientationParameters(...
        eop_file, ...
        initial_utc_time);

time_parameters.ut1_utc_sec                                   = ut1_utc_sec;
time_parameters.polar_motion_deg                              = polar_motion_deg;
time_parameters.longitude_earth_orientation_parameter_deg     = longitude_earth_orientation_parameter_deg;
time_parameters.obliquity_earth_orientation_parameter_deg     = obliquity_earth_orientation_parameter_deg;
time_parameters.length_of_day_earth_orientation_parameter_sec = length_of_day_earth_orientation_parameter_sec;

% Station positions
station_positions_itrf = transpose([
    -6143584, 1364250, 1033743;
    1907295, 6030810, -817119;
    2390310, -5564341, 1994578;
]);

%% Load the measurement
measurements                     = load("provided/LEO_DATA_Apparent-1.mat");
measurements.station_number      = measurements.LEO_DATA_Apparent(:, 1);
measurements.time                = measurements.LEO_DATA_Apparent(:, 2);
measurements.apparent_range      = measurements.LEO_DATA_Apparent(:, 3) * Units.kilometers;
measurements.apparent_range_rate = measurements.LEO_DATA_Apparent(:, 4) * Units.kilometers / Units.seconds;
measurements.LEO_DATA_Apparent   = [];

%% Load the trajectory
load("state_history.mat");

%% Iterate through the measurements and predict the associated measurements
num_measurements = numel(measurements.time);

predicted_measurements.time                 = measurements.time;
predicted_measurements.station_number       = measurements.station_number;
predicted_measurements.apparent_range       = zeros(num_measurements, 1);
predicted_measurements.apparent_range_rate  = zeros(num_measurements, 1);

integration_time = seconds(state_history.time - state_history.time(1));

for idx = 1:num_measurements
    % Find the state associated with the measurement
    measurement_time = measurements.time(idx);
    [~, state_idx] = min(abs(integration_time - measurement_time));
    assert(abs(integration_time(state_idx) - measurement_time) < 1e-6);

    % Extract the satellite's position and velocity
    satellite_position_at_rx_gcrf = state_history.satellite_position_gcrf(:, state_idx);
    satellite_velocity_at_rx_gcrf = state_history.satellite_velocity_gcrf(:, state_idx);

    station_position_at_rx_itrf = station_positions_itrf(:, measurements.station_number(idx));

    % Get the current GCRF/ITRF transform
    utc_time = state_history.time(state_idx);

    [C_itrf2gcrf, C_itrf2pef, C_pef2tod, C_tod2mod, C_mod2gcrf, R_dot, omega] = ...
        ComputeITRF2GCRF1976TransformParameters(...
            utc_time, ...
            time_parameters.ut1_utc_sec, ...
            time_parameters.polar_motion_deg, ...
            time_parameters.longitude_earth_orientation_parameter_deg, ...
            time_parameters.obliquity_earth_orientation_parameter_deg, ...
            time_parameters.length_of_day_earth_orientation_parameter_sec);

    % Transform the station position and velocity into the GCRF frame
    % Vallado, algorithm 24
    station_position_at_rx_gcrf = C_itrf2gcrf * station_position_at_rx_itrf;
    station_velocity_at_rx_gcrf = C_mod2gcrf * C_tod2mod * C_pef2tod * (cross([0;0;omega], C_itrf2pef * station_position_at_rx_itrf, 1));

    % Compute the apparent range and range-rate
    [apparent_range, apparent_range_rate, light_time] = ApparentRangeAndRangeRate( ...
        satellite_position_at_rx_gcrf, ...
        satellite_velocity_at_rx_gcrf, ...
        station_position_at_rx_gcrf, ...
        station_velocity_at_rx_gcrf, ...
        "threshold", 0.01 / Constants.SPEED_OF_LIGHT);

    predicted_measurements.apparent_range(idx) = apparent_range;
    predicted_measurements.apparent_range_rate(idx) = apparent_range_rate;
end

%% Now compute the difference
apparent_range_difference      = measurements.apparent_range      - predicted_measurements.apparent_range;
apparent_range_rate_difference = measurements.apparent_range_rate - predicted_measurements.apparent_range_rate;

rms_apparent_range_difference = rms(apparent_range_difference);
rms_apparent_range_rate_difference = rms(apparent_range_rate_difference);

figure('Position', [100, 100, 1250, 500]);
tlt = tiledlayout(1, 2);

nexttile(1);
scatter(measurements.time / 3600, apparent_range_difference / 1e3, 'ko', 'filled');
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Simulation Time [hours]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("Range Error [km]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("Apparent Range Errors (OMC)", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");

nexttile(2);
scatter(measurements.time / 3600, apparent_range_rate_difference / 1e3, 'ko', 'filled');
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Simulation Time [hours]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("Range Rate Error [km/s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("Apparent Range Rate Errors (OMC)", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");
