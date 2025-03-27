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

% Drag
drag_parameters.mass                      = 2000 .* Units.kilograms;
drag_parameters.reference_air_density     = 3.614e-13 .* Units.kilograms ./ Units.meters.^3;
drag_parameters.reference_radial_distance = gravity_parameters.earth_radius + 700000.0 .* Units.meters;
drag_parameters.decay_rate                = 88667.0 .* Units.meters;
drag_parameters.earth_rotation_rate       = 7.292115146706979e-5 .* Units.radians ./ Units.seconds;

% Satellite structure parameters
include_only_x_facet = true;
if include_only_x_facet
    structure_parameters.facet_descriptions = transpose([
        "+X";
    ]);
    
    structure_parameters.facet_orientations_body = transpose([
        +1, 0, 0;
    ]);
    
    structure_parameters.facet_areas = transpose([
        6;
    ]) .* (Units.meters).^2;
    
    structure_parameters.diffuse_reflection_coefficients = transpose([ ...
        0.04;
    ]);
    
    structure_parameters.specular_reflection_coefficients = transpose([ ...
        0.59;
    ]);
else
    structure_parameters.facet_descriptions = transpose([
        "+X", "-X", "+Y", "-Y", "+Z", "-Z", "Solar";
    ]);
    
    structure_parameters.facet_orientations_body = transpose([
        +1, 0, 0;
        -1, 0, 0;
        0, +1, 0;
        0, -1, 0;
        0, 0, +1;
        0, 0, -1;
        nan, nan, nan;
    ]);
    
    structure_parameters.facet_areas = transpose([
        6, 6, 8, 8, 12, 12, 15;
    ]) .* (Units.meters).^2;
    
    structure_parameters.diffuse_reflection_coefficients = transpose([ ...
        0.04, 0.04, 0.04 0.04, 0.80, 0.28, 0.04;
    ]);
    
    structure_parameters.specular_reflection_coefficients = transpose([ ...
        0.59, 0.59, 0.59, 0.59, 0.04, 0.18, 0.04;
    ]);
end


% Satellite state
initial_satellite_position_gcrf = [
    6990077.798814194; 
    1617465.311978378; 
    22679.810569245355;
] .* Units.meters;

initial_satellite_velocity_gcrf = [
    -1675.13972506056; 
    7273.72441330686; 
    252.688512916741;
] .* Units.meters / Units.seconds;

% Station positions
station_positions_itrf = transpose([
    -6143584, 1364250, 1033743;
    1907295, 6030810, -817119;
    2390310, -5564341, 1994578;
]);

% Coefficient of drag
initial_coefficient_of_drag = 2.0;

% Number of states in the filter
state_parameters.num_states = 7;

%% Perform the cubature transform
initial_state = [
    initial_satellite_position_gcrf;
    initial_satellite_velocity_gcrf;
    initial_coefficient_of_drag;
];

initial_error_covariance = blkdiag(1^2 * eye(3), 0.1^2 * eye(3), 0.001^2 * eye(1));

process_noise_covariance = zeros(7, 7);

cubature_transform_dynamics_function = @(state) CubatureTransformContinuousTimeDynamicsFunction(...
        initial_utc_time, ...
        state, ...
        state_parameters, ...
        time_parameters, ...
        gravity_parameters, ...
        drag_parameters, ...
        structure_parameters);

[d_dt_initial_state, ~, ~, A_cubature] = CubatureTransform( ...
    initial_state, ...
    initial_error_covariance, ...
    cubature_transform_dynamics_function, ...
    process_noise_covariance);

%% Compare against the provided A matrix
provided_A = load("provided/A_t0.mat");
provided_A = provided_A.A;

% The provided A matrix is in unit km and km / s. Make adjustments
provided_A(1:3, 7) = provided_A(1:3, 7) * 1e3;
provided_A(4:6, 7) = provided_A(4:6, 7) * 1e3;
provided_A(7, 1:3) = provided_A(7, 1:3) / 1e3;
provided_A(7, 4:6) = provided_A(7, 4:6) / 1e3;

%% Compare
A_absolute_error = A_cubature - provided_A;
A_relative_error = (A_cubature - provided_A) ./ provided_A;
