%% Preamble
clc; clear all; close all;

%% Provided Constants
earth_gravitational_parameter = 398600.4415 .* Units.kilometers.^3 ./ Units.seconds.^2;
earth_radius                  = 6378.1363 .* Units.kilometers;
sun_gravitational_parameter   = 132712440018 .* Units.kilometers.^3 ./ Units.seconds.^2;
moon_gravitational_parameter  = 4902.800066 .* Units.kilometers.^3 ./ Units.seconds.^2;
earth_eccentricity            = 0.081819221456;
earth_rotation_rate           = 7.292115146706979e-5 .* Units.radians ./ Units.seconds;

%% Provided
% Provided GCRF satellite position and velocity
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

% Provided station coordinates
station_positions_itrf = transpose([
    -6143584, 1364250, 1033743;
    1907295, 6030810, -817119;
    2390310, -5564341, 1994578;
]);

% Provided UTC date
utc_year    = 2018;
utc_month   = 02;
utc_day     = 1;
utc_hour    = 05;
utc_minute  = 00;
utc_seconds = 0;
utc_time    = datetime(utc_year, utc_month, utc_day, utc_hour, utc_minute, utc_seconds);

drag_model.mass                      = 2000 .* Units.kilograms;
drag_model.reference_air_density     = 3.614e-13 .* Units.kilograms ./ Units.meters.^3;
drag_model.reference_radial_distance = earth_radius + 700000.0 .* Units.meters;
drag_model.decay_rate                = 88667.0 .* Units.meters;

%% Parse the Earth orientation parameters
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

%% Convert time
converted_time = ConvertUTCTime(utc_time, ut1_utc_sec);

%% Compute the frame transform parameters
[C_itrf2gcrf, C_itrf2pef, C_pef2tod, C_tod2mod, C_mod2gcrf, R_dot, omega] = ...
    ComputeITRF2GCRF1976TransformParameters(...
        utc_time, ...
        ut1_utc_sec, ...
        polar_motion_deg, ...
        longitude_earth_orientation_parameter_deg, ...
        obliquity_earth_orientation_parameter_deg, ...
        length_of_day_earth_orientation_parameter_sec);
C_gcrf2itrf = transpose(C_itrf2gcrf);

%% Transform the station position and velocity into the GCRF frame
% Vallado, algorithm 24
station_positions_gcrf  = C_itrf2gcrf * station_positions_itrf;
station_velocities_gcrf = C_mod2gcrf * C_tod2mod * C_pef2tod * (cross(repmat([0;0;omega], 1, 3), C_itrf2pef * station_positions_itrf, 1));


%% Load the measurements
measurements                     = load("provided/LEO_DATA_Apparent-1.mat");
measurements.station_number      = measurements.LEO_DATA_Apparent(:, 1);
measurements.time                = measurements.LEO_DATA_Apparent(:, 2);
measurements.apparent_range      = measurements.LEO_DATA_Apparent(:, 3) * Units.kilometers;
measurements.apparent_range_rate = measurements.LEO_DATA_Apparent(:, 4) * Units.kilometers / Units.seconds;
measurements.LEO_DATA_Apparent   = [];

%% Load the provided Jacobian
load("provided/H_Tilde_t0.mat");

%% Predict the first measurement
satellite_position_at_rx_gcrf = satellite_position_gcrf;
satellite_velocity_at_rx_gcrf = satellite_velocity_gcrf;
station_position_at_rx_gcrf   = station_positions_gcrf(:, measurements.station_number(1));
station_velocity_at_rx_gcrf   = station_velocities_gcrf(:, measurements.station_number(1));

[apparent_range, apparent_range_rate, light_time] = ApparentRangeAndRangeRate( ...
        satellite_position_at_rx_gcrf, ...
        satellite_velocity_at_rx_gcrf, ...
        station_position_at_rx_gcrf, ...
        station_velocity_at_rx_gcrf, ...
        "threshold", 0.01 / Constants.SPEED_OF_LIGHT);

apparent_range_error      = apparent_range      - measurements.apparent_range(1);
apparent_range_rate_error = apparent_range_rate - measurements.apparent_range_rate(1);

disp(apparent_range_error);
disp(apparent_range_rate_error);

%% Compute the Jacobian of the first measurement
H = ApparentRangeAndRangeRateJacobian( ...
        satellite_position_at_rx_gcrf, ...
        satellite_velocity_at_rx_gcrf, ...
        station_position_at_rx_gcrf, ...
        station_velocity_at_rx_gcrf, ...
        light_time);

jacobian_error = H(1:2, 1:6) - H_TILDA(1:2, 1:6);
disp(jacobian_error);

%% Compute the Jacobian with the Cubature transform

function [z] = CubatureFunction(x, station_position_at_rx_gcrf, station_velocity_at_rx_gcrf)
    num_sample_points = size(x, 2);
    
    z = zeros(2, num_sample_points);
    
    for idx = 1:num_sample_points
        satellite_position_at_rx_gcrf = x(1:3, idx);
        satellite_velocity_at_rx_gcrf = x(4:6, idx);

        [apparent_range, apparent_range_rate, light_time] = ApparentRangeAndRangeRate( ...
                satellite_position_at_rx_gcrf, ...
                satellite_velocity_at_rx_gcrf, ...
                station_position_at_rx_gcrf, ...
                station_velocity_at_rx_gcrf, ...
                "threshold", 0.01 / Constants.SPEED_OF_LIGHT);

        z(1, idx) = apparent_range;
        z(2, idx) = apparent_range_rate;
    end
end

x = [satellite_position_at_rx_gcrf; satellite_velocity_at_rx_gcrf];
Pxx = blkdiag(1^2 * eye(3), 0.1^2 * eye(3));
h = @(x) CubatureFunction(x, station_position_at_rx_gcrf, station_velocity_at_rx_gcrf);
R = diag([1^2, 0.1^2]);
[~, ~, ~, H_cubature] = CubatureTransform(x, Pxx, h, R);


%% Compute the A matrix
satellite_position_itrf = C_gcrf2itrf * satellite_position_gcrf;

[gravitational_acceleration_jacobian_itrf] = ...
    EGM96GravitationalAccelerationJacobian(...
        satellite_position_itrf, ...
        earth_gravitational_parameter, ...
        earth_radius);

% Sandwich the Jacobian?
gravitational_acceleration_jacobian_gcrf = C_itrf2gcrf * gravitational_acceleration_jacobian_itrf * transpose(C_itrf2gcrf);


% coefficient_of_drag = 1.88;
coefficient_of_drag = 2;
effective_area      = 6 .* (Units.meters).^2;

drag_model.mass                      = 2000 .* Units.kilograms;
drag_model.reference_air_density     = 3.614e-13 .* Units.kilograms ./ Units.meters.^3;
drag_model.reference_radial_distance = earth_radius + 700000.0 .* Units.meters;
drag_model.decay_rate                = 88667.0 .* Units.meters;

[drag_position_gcrf_jacobian, drag_velocity_gcrf_jacobian, drag_coefficient_of_drag_jacobian] = ...
    SimpleDragJacobianGCRF(...
        satellite_position_gcrf, ...
        satellite_velocity_gcrf, ...
        coefficient_of_drag, ...
        effective_area, ...
        drag_model.mass, ...
        drag_model.reference_air_density, ...
        drag_model.reference_radial_distance, ...
        drag_model.decay_rate, ...
        earth_rotation_rate);

sun_position_tod = ComputeSunPositionTOD(converted_time.ut1, converted_time.tdb);
sun_position_gcrf = C_mod2gcrf * C_tod2mod * sun_position_tod;

moon_position_gcrf = ComputeMoonPositionGCRF( ...
        converted_time.tdb, ...
        earth_radius);

sun_position_jacobian_gcrf = ...
    SimpleTwoBodyOffsetGravitationalAccelerationJacobianGCRF(...
        (satellite_position_gcrf - sun_position_gcrf), ...
        sun_gravitational_parameter);

moon_position_jacobian_gcrf = ...
    SimpleTwoBodyOffsetGravitationalAccelerationJacobianGCRF(...
        (satellite_position_gcrf - moon_position_gcrf), ...
        moon_gravitational_parameter);

% d/dt p  = v;
% d/dt v  = gravity + drag + SRP + sun + moon
% d/dt CD = 0

A = zeros(7, 7);

A(1:3, 1:3) = zeros(3, 3);
A(1:3, 4:6) = eye(3);
A(1:3, 7)   = zeros(3, 1);

A(4:6, 1:3) = gravitational_acceleration_jacobian_gcrf + drag_position_gcrf_jacobian;
A(4:6, 4:6) = drag_velocity_gcrf_jacobian;
A(4:6, 7)   = drag_coefficient_of_drag_jacobian;

A(7, 1:3)   = zeros(1, 3);
A(7, 4:6)   = zeros(1, 3);
A(7, 7)     = 0;

%% Test
[drag_position_gcrf_jacobian1, drag_velocity_gcrf_jacobian1, drag_coefficient_of_drag_jacobian1] = ...
    SimpleDragJacobianGCRF(...
        satellite_position_gcrf, ...
        satellite_velocity_gcrf, ...
        1.88, ...
        effective_area, ...
        drag_model.mass, ...
        drag_model.reference_air_density, ...
        drag_model.reference_radial_distance, ...
        drag_model.decay_rate, ...
        earth_rotation_rate);

[drag_position_gcrf_jacobian2, drag_velocity_gcrf_jacobian2, drag_coefficient_of_drag_jacobian2] = ...
    SimpleDragJacobianGCRF(...
        satellite_position_gcrf, ...
        satellite_velocity_gcrf, ...
        2.0, ...
        effective_area, ...
        drag_model.mass, ...
        drag_model.reference_air_density, ...
        drag_model.reference_radial_distance, ...
        drag_model.decay_rate, ...
        earth_rotation_rate);

%% Compare against the provided A matrix
provided_A = load("provided/A_t0.mat");
provided_A = provided_A.A;

% The provided A matrix is in unit km and km / s. Make adjustments
provided_A(1:3, 7) = provided_A(1:3, 7) * 1e3;
provided_A(4:6, 7) = provided_A(4:6, 7) * 1e3;
provided_A(7, 1:3) = provided_A(7, 1:3) / 1e3;
provided_A(7, 4:6) = provided_A(7, 4:6) / 1e3;


%% Transform position into GCRF frame
% position_pef  = C_itrf2pef * position_itrf;
% position_tod  = C_pef2tod  * position_pef;
% position_mod  = C_tod2mod  * position_tod;
% position_gcrf = C_mod2gcrf * position_mod;
% 
% velocity_pef  = C_itrf2pef * velocity_itrf;
% velocity_tod  = C_pef2tod  * velocity_pef + R_dot * position_pef;
% velocity_mod  = C_tod2mod  * velocity_tod;
% velocity_gcrf = C_mod2gcrf * velocity_mod;


%% Dynamics function
time_parameters.ut1_utc_sec                                   = ut1_utc_sec;
time_parameters.polar_motion_deg                              = polar_motion_deg;
time_parameters.longitude_earth_orientation_parameter_deg     = longitude_earth_orientation_parameter_deg;
time_parameters.obliquity_earth_orientation_parameter_deg     = obliquity_earth_orientation_parameter_deg;
time_parameters.length_of_day_earth_orientation_parameter_sec = length_of_day_earth_orientation_parameter_sec;

gravity_parameters.earth_gravitational_parameter = earth_gravitational_parameter;
gravity_parameters.earth_radius                  = earth_radius;
gravity_parameters.spherical_harmonic_degree     = 20;
gravity_parameters.moon_gravitational_parameter  = moon_gravitational_parameter;
gravity_parameters.sun_gravitational_parameter   = sun_gravitational_parameter;

drag_parameters.effective_area            = 6 .* (Units.meters).^2; % TODO
drag_parameters.mass                      = drag_model.mass;
drag_parameters.reference_air_density     = drag_model.reference_air_density;
drag_parameters.reference_radial_distance = drag_model.reference_radial_distance;
drag_parameters.decay_rate                = drag_model.decay_rate;
drag_parameters.earth_rotation_rate       = earth_rotation_rate;

[d_dt_satellite_position_gcrf, d_dt_satellite_velocity_gcrf, d_dt_coefficient_of_drag] = ...
    DynamicsFunction(...
        utc_time, ...
        satellite_position_gcrf, ...
        satellite_velocity_gcrf, ...
        coefficient_of_drag, ...
        time_parameters, ...
        gravity_parameters, ...
        drag_parameters);

function [dx_dt] = DynamicsCubatureFunction(...
    x, ...
    utc_time, ...
    time_parameters, ...
    gravity_parameters, ...
    drag_parameters)

    % Unpack the state samples
    satellite_position_gcrf = x(1:3, :);
    satellite_velocity_gcrf = x(4:6, :);
    coefficient_of_drag     = x(7, :);

    % Push the state samples through the dynamics function
    [d_dt_satellite_position_gcrf, d_dt_satellite_velocity_gcrf, d_dt_coefficient_of_drag] = ...
        DynamicsFunction(...
            utc_time, ...
            satellite_position_gcrf, ...
            satellite_velocity_gcrf, ...
            coefficient_of_drag, ...
            time_parameters, ...
            gravity_parameters, ...
            drag_parameters);

    % Pack the state samples
    dx_dt = [d_dt_satellite_position_gcrf; d_dt_satellite_velocity_gcrf; reshape(d_dt_coefficient_of_drag, 1, [])];
end

x = [
    satellite_position_gcrf;
    satellite_velocity_gcrf;
    coefficient_of_drag;
];

Pxx = blkdiag(1^2 * eye(3), 0.1^2 * eye(3), 0.001^2 * eye(1));

Q   = zeros(7, 7);

f = @(x) DynamicsCubatureFunction(...
        x, ...
        utc_time, ...
        time_parameters, ...
        gravity_parameters, ...
        drag_parameters);

[dx_dt, ~, ~, A_cubature] = CubatureTransform(x, Pxx, f, Q);