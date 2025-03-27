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
include_only_x_facet = false;
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

%% Load truth data
% Truth position, velocity time
truth_data                               = load("provided/lighttime_truth-1.mat");
truth_data.rx_time                       = truth_data.lighttime_truth(:, 7);
truth_data.satellite_position_at_tx_gcrf = transpose(truth_data.lighttime_truth(:, 1:3)) * Units.kilometers;
truth_data.satellite_velocity_at_tx_gcrf = transpose(truth_data.lighttime_truth(:, 4:6)) * Units.kilometers / Units.seconds;
truth_data.lighttime_truth               = [];

%% Set up state
% Satellite state
initial_satellite_position_gcrf = truth_data.satellite_position_at_tx_gcrf(:, 1);
initial_satellite_velocity_gcrf = truth_data.satellite_velocity_at_tx_gcrf(:, 1);

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

%% Propagate the state and error covariance with the cubature transform
% Set up cubature parameters
% TODO: Configure initial error covariance and process noise covariance
initial_error_covariance = blkdiag((1)^2 * eye(3), 0.01^2 * eye(3), 0.1^2 * eye(1));
process_noise_covariance = zeros(state_parameters.num_states, state_parameters.num_states);

% Set up time
% TODO: Technically this is the satellite position at the time of RX, but since we're not using any measurements, we'll
% just use it as normal time.
time_history   = truth_data.rx_time + DateTimeToIntegrationTime(initial_utc_time);
time_history   = time_history(1:7);
initial_time   = time_history(1);
num_time_steps = numel(time_history);

% Set up state
initial_state = [
    initial_satellite_position_gcrf;
    initial_satellite_velocity_gcrf;
    initial_coefficient_of_drag;
];
state_history       = zeros(state_parameters.num_states, num_time_steps);
state_history(:, 1) = initial_state;

% Set up the error covariance
error_covariance_history          = zeros(state_parameters.num_states, state_parameters.num_states, num_time_steps);
error_covariance_history(:, :, 1) = initial_error_covariance;

tic;
for idx = 1:(num_time_steps - 1)
    % Get the current state
    current_integration_time = time_history(idx);
    current_state            = state_history(:, idx);
    current_error_covariance = error_covariance_history(:, :, idx);
    
    % Get the next time
    propagated_integration_time = time_history(idx + 1);

    % Configure the integration parameters
    integration_parameters.integration_start_time = current_integration_time;
    integration_parameters.integration_end_time   = propagated_integration_time;
    integration_parameters.ode45_options          = odeset('reltol', 3e-14, 'abstol', 1e-16);

    % This is what the Cubature transform uses. It's a wrapper for a function that integrates the dynamics function with
    % ode45.
    cubature_transform_discrete_time_dynamics_function = ...
        @(state) CubatureTransformDiscreteTimeDynamicsFunction(...
            state, ...
            state_parameters, ...
            time_parameters, ...
            gravity_parameters, ...
            drag_parameters, ...
            structure_parameters, ...
            integration_parameters);

    [propagated_state, propagated_error_covariance, ~, ~] = CubatureTransform( ...
        current_state, ...
        current_error_covariance, ...
        cubature_transform_discrete_time_dynamics_function, ...
        process_noise_covariance);

    state_history(:, idx + 1)               = propagated_state;
    error_covariance_history(:, :, idx + 1) = propagated_error_covariance;
end
toc;


%% Check the errors
position_error = state_history(1:3, :) - truth_data.satellite_position_at_tx_gcrf(:, 1:7);
velocity_error = state_history(4:6, :) - truth_data.satellite_velocity_at_tx_gcrf(:, 1:7);

% Gravitation only
% -17.5715594133362        -0.0265213392740407
%  25.7219295697287        -0.0178243468189976
%   1.0133837688627        -0.000343455862406472

% Gravitation and X-drag
% -17.5711009781808       -0.0265187283148407
%  25.7208255361766       -0.0178303818156564
%  1.01333990183775     -0.000343696112338421

% Gravitation and full drag
% -17.5702359694988        -0.026512818856645
% 25.7183116525412       -0.0178441033967829
% 1.01324136002222     -0.000344241720938498

% Gravitation and moon
% -17.5565704852343       -0.0264675114303827
%  25.6634719059803       -0.0181394482015094
% 0.993777019524714     -0.000441341076225399

% Gravitation and sun
% -17.5747442077845       -0.0265509055752773
%  25.7000829153694       -0.0179397047340899
%  1.00668761893758     -0.000375149295479105

% Gravitation and sun and moon
% -17.5598286110908       -0.0264970690004702
%  25.6415759730153       -0.0182549002047381
% 0.987079615922994     -0.000473037801384635

% Gravitation and sun and moon and full drag
% -17.5585182383657       -0.0264885513388435
%  25.6379497996531       -0.0182746723576201
% 0.986937001682236      -0.00047382418097186

% Gravitaiton sun moon and full drag
% -17.5819293344393       -0.0265370924234958
%  25.7986898240633       -0.0174135042261696
%  1.03954613261158     -0.000214658184347627
