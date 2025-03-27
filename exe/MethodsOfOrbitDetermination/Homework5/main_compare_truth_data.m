%% Preamble
clc; clear all; close all;

%% Provided Parameters
% Gravity
gravity_parameters.earth_gravitational_parameter = 398600.4415 .* Units.kilometers.^3 ./ Units.seconds.^2;
gravity_parameters.earth_radius                  = 6378.1363 .* Units.kilometers;
gravity_parameters.sun_gravitational_parameter   = 132712440018 .* Units.kilometers.^3 ./ Units.seconds.^2;
gravity_parameters.moon_gravitational_parameter  = 4902.800066 .* Units.kilometers.^3 ./ Units.seconds.^2;
gravity_parameters.spherical_harmonic_degree     = 20;

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

%% Compare truth data 
% Truth position, velocity time
truth_data                               = load("provided/lighttime_truth-1.mat");
truth_data.rx_time                       = truth_data.lighttime_truth(:, 7);
truth_data.satellite_position_at_tx_gcrf = transpose(truth_data.lighttime_truth(:, 1:3)) * Units.kilometers;
truth_data.satellite_velocity_at_tx_gcrf = transpose(truth_data.lighttime_truth(:, 4:6)) * Units.kilometers / Units.seconds;
truth_data.lighttime_truth               = [];

% Provided measurements
measurements                     = load("provided/LEO_DATA_Apparent-1.mat");
measurements.station_number      = measurements.LEO_DATA_Apparent(:, 1);
measurements.rx_time             = measurements.LEO_DATA_Apparent(:, 2);
measurements.apparent_range      = measurements.LEO_DATA_Apparent(:, 3) * Units.kilometers;
measurements.apparent_range_rate = measurements.LEO_DATA_Apparent(:, 4) * Units.kilometers / Units.seconds;
measurements.LEO_DATA_Apparent   = [];

num_measurements = numel(truth_data.rx_time);

true_apparent_range      = zeros(num_measurements, 1);
true_apparent_range_rate = zeros(num_measurements, 1);
true_light_time          = zeros(num_measurements, 1);

for idx = 1:num_measurements

    % Compute the time of reception
    dt_seconds = measurements.rx_time(idx);
    converted_time = ConvertUTCTime( ...
        initial_utc_time + seconds(dt_seconds), ...
        time_parameters.ut1_utc_sec);
    
    % Compute the ITRF to GCRF transform
    [C_itrf2gcrf, C_itrf2pef, C_pef2tod, C_tod2mod, C_mod2gcrf, R_dot, omega] = ...
        ComputeITRF2GCRF1976TransformParameters(...
            converted_time.utc, ...
            time_parameters.ut1_utc_sec, ...
            time_parameters.polar_motion_deg, ...
            time_parameters.longitude_earth_orientation_parameter_deg, ...
            time_parameters.obliquity_earth_orientation_parameter_deg, ...
            time_parameters.length_of_day_earth_orientation_parameter_sec);

    % Transform the station position into the GCRF frame
    station_positions_at_rx_gcrf  = C_itrf2gcrf * station_positions_itrf;
    station_velocities_at_rx_gcrf = C_mod2gcrf * C_tod2mod * C_pef2tod * (cross(repmat([0;0;omega], 1, 3), C_itrf2pef * station_positions_itrf, 1));

    station_idx                 = measurements.station_number(idx);
    station_position_at_rx_gcrf = station_positions_at_rx_gcrf(:, station_idx);
    station_velocity_at_rx_gcrf = station_velocities_at_rx_gcrf(:, station_idx);

    true_satellite_position_at_tx_gcrf = truth_data.satellite_position_at_tx_gcrf(:, idx);
    true_satellite_velocity_at_tx_gcrf = truth_data.satellite_velocity_at_tx_gcrf(:, idx);

    relative_position  = station_position_at_rx_gcrf - true_satellite_position_at_tx_gcrf;
    relative_velocity  = station_velocity_at_rx_gcrf - true_satellite_velocity_at_tx_gcrf;
    relative_direction = relative_position ./ norm(relative_position);

    apparent_range      = vecnorm(relative_position, 2, 1);
    apparent_range_rate = dot(relative_velocity, relative_direction);

    true_apparent_range(idx)      = apparent_range;
    true_apparent_range_rate(idx) = apparent_range_rate;
end

measured_apparent_range = measurements.apparent_range;
measured_apparent_range_rate = measurements.apparent_range_rate;

measured_apparent_range_error      = measured_apparent_range - true_apparent_range;
measured_apparent_range_rate_error = measured_apparent_range_rate - true_apparent_range_rate;

