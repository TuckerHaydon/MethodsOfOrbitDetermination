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
utc_year    = 2018;
utc_month   = 02;
utc_day     = 1;
utc_hour    = 05;
utc_minute  = 00;
utc_seconds = 0;
utc_time    = datetime(...
    utc_year, ...
    utc_month, ...
    utc_day, ...
    utc_hour, ...
    utc_minute, ...
    utc_seconds);

eop_file = "finals.all.iau1980.txt";
[...      
        ut1_utc_sec, ...
        polar_motion_deg, ...
        length_of_day_earth_orientation_parameter_sec, ...
        longitude_earth_orientation_parameter_deg, ...
        obliquity_earth_orientation_parameter_deg] = ...
    ParseEarthOrientationParameters(...
        eop_file, ...
        utc_time);

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
satellite_position_gcrf = [
    6990077.798814194; 
    1617465.311978378; 
    22679.810569245355;
] .* Units.meters;

satellite_velocity_gcrf = [
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
coefficient_of_drag = 2.0;

% Number of states in the filter
state_parameters.num_states = 7;

%% Preprocessing
% (1) Time
converted_time = ConvertUTCTime( ...
    utc_time, ...
    time_parameters.ut1_utc_sec);

% (2) Coordinate frames
[C_itrf2gcrf, C_itrf2pef, C_pef2tod, C_tod2mod, C_mod2gcrf, R_dot, omega] = ...
    ComputeITRF2GCRF1976TransformParameters(...
        utc_time, ...
        time_parameters.ut1_utc_sec, ...
        time_parameters.polar_motion_deg, ...
        time_parameters.longitude_earth_orientation_parameter_deg, ...
        time_parameters.obliquity_earth_orientation_parameter_deg, ...
        time_parameters.length_of_day_earth_orientation_parameter_sec);

C_gcrf2itrf = transpose(C_itrf2gcrf);

satellite_position_itrf = C_gcrf2itrf * satellite_position_gcrf;

% (3) Sun/Moon Position
sun_position_tod = ComputeSunPositionTOD( ...
    converted_time.ut1, ...
    converted_time.tdb);
sun_position_gcrf = C_mod2gcrf * C_tod2mod * sun_position_tod;

moon_position_tod = ComputeMoonPositionTOD( ...
    converted_time.tdb, ...
    gravity_parameters.earth_radius);
moon_position_gcrf = C_mod2gcrf * C_tod2mod * moon_position_tod;

% (4) Effective drag area
relative_sun_position_gcrf  = sun_position_gcrf - satellite_position_gcrf;
relative_sun_direction_gcrf = relative_sun_position_gcrf ./ vecnorm(relative_sun_position_gcrf, 2, 1);
velocity_direction_gcrf     = satellite_velocity_gcrf ./ vecnorm(satellite_velocity_gcrf, 2, 1);

% Transform from Transverse-Normal-Radial (body XYZ) to GCRF.
C_body2gcrf = AttitudeTNRToGCRF( ...
        satellite_position_gcrf, ...
        satellite_velocity_gcrf);

% Transform facet orientations from body frame to GCRF frame.
facet_orientations_gcrf = pagemtimes(...
    C_body2gcrf, ...
    structure_parameters.facet_orientations_body);

% Add solar panel facet orientation.
solar_panel_facet_mask = strcmp(structure_parameters.facet_descriptions, "Solar");
if any(solar_panel_facet_mask)
    facet_orientations_gcrf(:, solar_panel_facet_mask, :) = relative_sun_direction_gcrf;
end

facet_areas = structure_parameters.facet_areas;

% Ignore the fact that the solar panel is backwards in the drag computation
ignore_facet_backwards = false(size(facet_areas));
if any(solar_panel_facet_mask)
    ignore_facet_backwards(solar_panel_facet_mask, :) = true;
end

% Compute the effective area subject to drag
effective_drag_area = ComputeEffectiveDragArea( ...
        velocity_direction_gcrf, ...
        facet_orientations_gcrf, ...
        facet_areas, ...
        ignore_facet_backwards);

%% Evluate the initial A matrix with the Jacobians
earth_gravitational_acceleration_position_jacobian_itrf = ...
    EGM96GravitationalAccelerationJacobian(...
        satellite_position_itrf, ...
        gravity_parameters.earth_gravitational_parameter, ...
        gravity_parameters.earth_radius);

earth_gravitational_acceleration_position_jacobian_gcrf = ...
    C_itrf2gcrf * earth_gravitational_acceleration_position_jacobian_itrf * transpose(C_itrf2gcrf);

[drag_position_gcrf_jacobian, drag_velocity_gcrf_jacobian, drag_coefficient_of_drag_jacobian] = ...
    SimpleDragJacobianGCRF(...
        satellite_position_gcrf, ...
        satellite_velocity_gcrf, ...
        coefficient_of_drag, ...
        effective_drag_area, ...
        drag_parameters.mass, ...
        drag_parameters.reference_air_density, ...
        drag_parameters.reference_radial_distance, ...
        drag_parameters.decay_rate, ...
        drag_parameters.earth_rotation_rate);

moon_gravitational_acceleration_position_jacobian_gcrf = ...
    ThirdBodyPerturbingAccelerationJacobianGCRF(...
        satellite_position_gcrf, ...
        moon_position_gcrf, ...
        gravity_parameters.moon_gravitational_parameter);

sun_gravitational_position_jacobian_gcrf = ...
    ThirdBodyPerturbingAccelerationJacobianGCRF(...
        satellite_position_gcrf, ...
        sun_position_gcrf, ...
        gravity_parameters.sun_gravitational_parameter);

analytical_A = zeros(7, 7);

analytical_A(1:3, 1:3) = zeros(3, 3);
analytical_A(1:3, 4:6) = eye(3);
analytical_A(1:3, 7)   = zeros(3, 1);

analytical_A(4:6, 1:3) = ...
    earth_gravitational_acceleration_position_jacobian_gcrf + ...
    drag_position_gcrf_jacobian + ...
    moon_gravitational_acceleration_position_jacobian_gcrf + ...
    sun_gravitational_position_jacobian_gcrf;
analytical_A(4:6, 4:6) = drag_velocity_gcrf_jacobian;
analytical_A(4:6, 7)   = drag_coefficient_of_drag_jacobian;

analytical_A(7, 1:3) = zeros(1, 3);
analytical_A(7, 4:6) = zeros(1, 3);
analytical_A(7, 7)   = 0;

%% Evaluate the initial A matrix with the cubature transform
initial_state = [
    satellite_position_gcrf;
    satellite_velocity_gcrf;
    coefficient_of_drag;
];

initial_error_covariance = blkdiag(1^2 * eye(3), 0.1^2 * eye(3), 0.001^2 * eye(1));

process_noise_covariance = zeros(7, 7);

cubature_transform_dynamics_function = @(state) CubatureTransformContinuousTimeDynamicsFunction(...
        utc_time, ...
        state, ...
        state_parameters, ...
        time_parameters, ...
        gravity_parameters, ...
        drag_parameters, ...
        structure_parameters);

[d_dt_initial_state, ~, ~, cubature_A] = CubatureTransform( ...
    initial_state, ...
    initial_error_covariance, ...
    cubature_transform_dynamics_function, ...
    process_noise_covariance);

%% Load the provided A matrix
provided_A = load("provided/A_t0.mat");
provided_A = provided_A.A;

% The provided A matrix is in unit km and km / s. Make adjustments
provided_A(1:3, 7) = provided_A(1:3, 7) * 1e3;
provided_A(4:6, 7) = provided_A(4:6, 7) * 1e3;
provided_A(7, 1:3) = provided_A(7, 1:3) / 1e3;
provided_A(7, 4:6) = provided_A(7, 4:6) / 1e3;

%% Compare
analytical_A_absolute_error = analytical_A - provided_A;
analytical_A_relative_error = (analytical_A - provided_A) ./ provided_A;

cubature_A_absolute_error = cubature_A - provided_A;
cubature_A_relative_error = (cubature_A - provided_A) ./ provided_A;

%% Make plot
figure('Position', [100, 100, 600, 400]);
histogram(reshape(log10(abs(analytical_A_relative_error(4:7, 1:7))),[],1), 10);
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Log10 Relative Error", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("Count", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("Relative Error Magnitude of Jacobian Dynamics", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");

%% Load the measurement
measurements                     = load("provided/LEO_DATA_Apparent-1.mat");
measurements.station_number      = measurements.LEO_DATA_Apparent(:, 1);
measurements.time                = measurements.LEO_DATA_Apparent(:, 2);
measurements.apparent_range      = measurements.LEO_DATA_Apparent(:, 3) * Units.kilometers;
measurements.apparent_range_rate = measurements.LEO_DATA_Apparent(:, 4) * Units.kilometers / Units.seconds;
measurements.LEO_DATA_Apparent   = [];

%% Compute the measurement Jacobian
% Transform the station position and velocity into the GCRF frame
% Vallado, algorithm 24
station_positions_gcrf  = C_itrf2gcrf * station_positions_itrf;
station_velocities_gcrf = C_mod2gcrf * C_tod2mod * C_pef2tod * (cross(repmat([0;0;omega], 1, 3), C_itrf2pef * station_positions_itrf, 1));

[apparent_range, apparent_range_rate, light_time] = ApparentRangeAndRangeRate( ...
        satellite_position_gcrf, ...
        satellite_velocity_gcrf, ...
        station_positions_gcrf(:, measurements.station_number(1)), ...
        station_velocities_gcrf(:, measurements.station_number(1)), ...
        "threshold", 0.01 / Constants.SPEED_OF_LIGHT);

H = ApparentRangeAndRangeRateJacobian( ...
        satellite_position_gcrf, ...
        satellite_velocity_gcrf, ...
        station_positions_gcrf(:, measurements.station_number(1)), ...
        station_velocities_gcrf(:, measurements.station_number(1)), ...
        light_time);

%% Load the provided H matrix
provided_H = load("provided/H_Tilde_t0.mat");
provided_H = provided_H.H_TILDA;
provided_H = provided_H(:, 1:6);

H_absolute_error = H - provided_H;
H_relative_error = (H - provided_H) ./ provided_H;