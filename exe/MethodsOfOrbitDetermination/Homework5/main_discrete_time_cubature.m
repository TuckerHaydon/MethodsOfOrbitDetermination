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

%% Set up ODE 45 and Cubature transform
function [d_dt_state] = ODE45DynamicsFunction(...
        integration_time, ...
        state, ...
        state_parameters, ...
        time_parameters, ...
        gravity_parameters, ...
        drag_parameters, ...
        structure_parameters)
    
    % Reshape the state so its a [num_states, N] set of appended state vectors.
    state = reshape(state, state_parameters.num_states, []);

    % Unpack the state samples
    satellite_position_gcrf = state(1:3, :);
    satellite_velocity_gcrf = state(4:6, :);
    coefficient_of_drag     = state(7, :);

    % Convert integration time to UTC date time
    utc_time = IntegrationTimeToDateTime(integration_time);

    % Push the state samples through the dynamics function
    [d_dt_satellite_position_gcrf, d_dt_satellite_velocity_gcrf, d_dt_coefficient_of_drag] = ...
        DynamicsFunction(...
            utc_time, ...
            satellite_position_gcrf, ...
            satellite_velocity_gcrf, ...
            coefficient_of_drag, ...
            time_parameters, ...
            gravity_parameters, ...
            drag_parameters, ...
            structure_parameters);

    % Pack the state samples
    d_dt_state = [
        d_dt_satellite_position_gcrf; 
        d_dt_satellite_velocity_gcrf; 
        reshape(d_dt_coefficient_of_drag, 1, []);
    ];

    % Reshape back into a column vector
    d_dt_state = d_dt_state(:);
end

function [propagated_state] = CubatureTransformDiscreteTimeDynamicsFunction(...
        state, ...
        state_parameters, ...
        time_parameters, ...
        gravity_parameters, ...
        drag_parameters, ...
        structure_parameters, ...
        integration_parameters)
    arguments(Input)
        state(:, :) double
        state_parameters(1, 1) struct
        time_parameters(1, 1) struct
        gravity_parameters(1, 1) struct
        drag_parameters(1, 1) struct
        structure_parameters(1, 1) struct
        integration_parameters(1, 1) struct
    end

    arguments(Output)
        propagated_state(:, :) double
    end

    % State should be size [num_states, N]
    assert(size(state, 1) == state_parameters.num_states);

    % Convert time to modified Julian date in seconds.
    start_time = integration_parameters.integration_start_time;
    end_time   = integration_parameters.integration_end_time;
    time_span  = [start_time, end_time];

    % Set up the dynamics function
    ode_45_dynamics_function = ...
        @(t, y) ODE45DynamicsFunction(...
            t, ...
            y, ...
            state_parameters, ...
            time_parameters, ...
            gravity_parameters, ...
            drag_parameters, ...
            structure_parameters);

    % ODE 45 requires the state to be a vector.
    state = state(:);

    % Integrate the state forward in time.
    [~, propagated_state_history] = ode45( ...
        ode_45_dynamics_function, ...
        time_span, ...
        state, ...
        integration_parameters.ode45_options);

    % Each stacked state is a row, and the final row is the final state. Extract it and reshape it.
    propagated_state = reshape(transpose(propagated_state_history(end, :)), state_parameters.num_states, []);

    % State should be size [num_states, N]
    assert(size(propagated_state, 1) == state_parameters.num_states);
end

%% Propagate the state and error covariance with the cubature transform
% Set up cubature parameters
% TODO: Configure initial error covariance and process noise covariance
initial_error_covariance = blkdiag((2e3)^2 * eye(3), 20^2 * eye(3), 0.1^2 * eye(1));
process_noise_covariance = zeros(state_parameters.num_states, state_parameters.num_states);

% Set up time
initial_time   = DateTimeToIntegrationTime(initial_utc_time);
final_time     = initial_time + 600 .* Units.seconds;
time_step      = 60 .* Units.seconds;
time_history   = reshape(colon(initial_time, time_step, final_time), [], 1);
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


