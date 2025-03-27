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

%% Integrate the state
initial_utc_integration_time = initial_utc_time + seconds(0);
final_utc_integration_time   = initial_utc_time + seconds(21600);

time_span = colon(...
    DateTimeToIntegrationTime(initial_utc_integration_time), ...
    60, ...
    DateTimeToIntegrationTime(final_utc_integration_time));

% ode45_options = odeset('reltol', 3e-14, 'abstol', 1e-16);
ode45_options = odeset('reltol', 3e-11, 'abstol', 1e-13);

initial_stm = eye(7, 7);

state_vector = [
    initial_satellite_position_gcrf;
    initial_satellite_velocity_gcrf;
    initial_coefficient_of_drag;
    initial_stm(:);
];

function [y_dot] = ODE45StateAndSTMDynamicsFunction(...
        t, ...
        y, ...
        time_parameters, ...
        gravity_parameters, ...
        drag_parameters, ...
        structure_parameters)
    
    utc_time = IntegrationTimeToDateTime(t);

    satellite_position_gcrf = y(1:3);
    satellite_velocity_gcrf = y(4:6);
    coefficient_of_drag     = y(7);
    state_transition_matrix = reshape(y(8:end), 7, 7);


    [d_dt_satellite_position_gcrf, d_dt_satellite_velocity_gcrf, d_dt_coefficient_of_drag, d_dt_state_transition_matrix] = ...
        StateAndSTMDynamicsFunction(...
            utc_time, ...
            satellite_position_gcrf, ...
            satellite_velocity_gcrf, ...
            coefficient_of_drag, ...
            state_transition_matrix, ...
            time_parameters, ...
            gravity_parameters, ...
            drag_parameters, ...
            structure_parameters);

    y_dot = [
        d_dt_satellite_position_gcrf;
        d_dt_satellite_velocity_gcrf;
        d_dt_coefficient_of_drag;
        d_dt_state_transition_matrix(:);
    ];
end

% Set up the dynamics function
ode_45_dynamics_function = ...
    @(t, y) ODE45StateAndSTMDynamicsFunction(...
        t, ...
        y, ...
        time_parameters, ...
        gravity_parameters, ...
        drag_parameters, ...
        structure_parameters);

% Integrate the state forward in time.
% profile on
tic;
[~, propagated_state_history] = ode45( ...
    ode_45_dynamics_function, ...
    time_span, ...
    state_vector, ...
    ode45_options);
toc;
% profile viewer

final_state = propagated_state_history(end, 1:7);
final_stm   = reshape(propagated_state_history(end, 8:end), 7, 7);

%% Load the provided STM
provided_stm = load("provided/Phi_21600_0(1).mat");
provided_stm = provided_stm.PHI_t_120;

% Make km to meters adjustment
provided_stm(1:6, 7) = provided_stm(1:6, 7) * 1e3;
provided_stm(7, 1:6) = provided_stm(1:6, 7) / 1e3;

%% Compute the relative error of the block 1
relative_stm_error = (final_stm - provided_stm) ./ provided_stm;
relative_position_stm_error = relative_stm_error(1:3, 1:3);

figure('Position', [100, 100, 600, 400]);
histogram(reshape(log10(abs(relative_stm_error(1:6, :))),[],1), 10);
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Log10 Relative Error", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("Count", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("Relative Error Magnitude of State Transition Matrix", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");

%% Plot the trajectory
figure('Position', [100, 100, 900, 900]);
plot3(...
    propagated_state_history(:, 1) / 1e3, ...
    propagated_state_history(:, 2) / 1e3, ...
    propagated_state_history(:, 3) / 1e3, ...
    'k-', "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("GCRF X [km]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("GCRF Y [km]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
zlabel("GCRF Z [km]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("Integrated Orbit Over 6 Hours", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");

%% Load the truth trajectory
% Truth position, velocity time
truth_data                               = load("provided/lighttime_truth-1.mat");
truth_data.rx_time                       = truth_data.lighttime_truth(:, 7);
truth_data.satellite_position_at_tx_gcrf = transpose(truth_data.lighttime_truth(:, 1:3)) * Units.kilometers;
truth_data.satellite_velocity_at_tx_gcrf = transpose(truth_data.lighttime_truth(:, 4:6)) * Units.kilometers / Units.seconds;
truth_data.lighttime_truth               = [];

%% Compute and plot the trajectory error
time = time_span - time_span(1);
state_mask = false(size(propagated_state_history, 1), 1);
for idx = 1:numel(state_mask)
    if any(abs(truth_data.rx_time - time(idx)) < 1e-6)
        state_mask(idx) = true;
    end
end
assert(sum(state_mask) == numel(truth_data.rx_time));

satellite_position_error_gcrf = transpose(propagated_state_history(state_mask, 1:3)) - truth_data.satellite_position_at_tx_gcrf;
satellite_velocity_error_gcrf = transpose(propagated_state_history(state_mask, 4:6)) - truth_data.satellite_velocity_at_tx_gcrf;

figure('Position', [100, 100, 1250, 500]);
tlt = tiledlayout(1, 2);

nexttile(1);
h = plot(truth_data.rx_time / 3600, satellite_position_error_gcrf / 1e3, '-', 'LineWidth', 3.0);
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Simulation Time [hours]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("GCRF Position Error [km]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("GCRF Position Error", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");
legend([h(1), h(2), h(3)], "GCRF X", "GCRF Y", "GCRF Z", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");

nexttile(2);
h = plot(truth_data.rx_time / 3600, satellite_velocity_error_gcrf / 1e3, '-', 'LineWidth', 3.0);
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Simulation Time [hours]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("GCRF Velocity Error [km/s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("GCRF Velocity Error", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");
legend([h(1), h(2), h(3)], "GCRF X", "GCRF Y", "GCRF Z", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");


%% Save the integrated state
state_history.time = IntegrationTimeToDateTime(time_span);
state_history.satellite_position_gcrf = transpose(propagated_state_history(:, 1:3));
state_history.satellite_velocity_gcrf = transpose(propagated_state_history(:, 4:6));
state_history.coefficient_of_drag     = propagated_state_history(:, 7);

save("state_history.mat", "state_history");
