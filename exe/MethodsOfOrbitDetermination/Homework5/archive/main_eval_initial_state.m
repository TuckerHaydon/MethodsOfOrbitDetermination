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
utc_time    = datetime(utc_year, utc_month, utc_day, utc_hour, utc_minute, utc_seconds);

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


%% Eval initial dynamics
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

initial_satellite_position_itrf = C_gcrf2itrf * initial_satellite_position_gcrf;

%% Dynamics
% (1) Earth Gravity
earth_gravitational_acceleration_itrf = EGM96GravitationalAcceleration(...
    initial_satellite_position_itrf, ...
    gravity_parameters.earth_gravitational_parameter, ...
    gravity_parameters.earth_radius, ...
    gravity_parameters.spherical_harmonic_degree);

earth_gravitational_acceleration_gcrf = C_itrf2gcrf * earth_gravitational_acceleration_itrf;

% (2) Drag
% (2a) Compute the position of the Sun and Moon
sun_position_tod = ComputeSunPositionTOD( ...
    converted_time.ut1, ...
    converted_time.tdb);
sun_position_gcrf = C_mod2gcrf * C_tod2mod * sun_position_tod;

moon_position_tod = ComputeMoonPositionTOD( ...
    converted_time.tdb, ...
    gravity_parameters.earth_radius);
moon_position_gcrf = C_mod2gcrf * C_tod2mod * moon_position_tod;

% (2b) Compute the Sun direction (the solar panel is pointed towards it)
relative_sun_position_gcrf      = sun_position_gcrf - initial_satellite_position_gcrf;
relative_sun_direction_gcrf     = relative_sun_position_gcrf ./ norm(relative_sun_position_gcrf);
initial_velocity_direction_gcrf = initial_satellite_velocity_gcrf ./ norm(initial_satellite_velocity_gcrf);

% (2c) Compute the facet orientations in the GCRF frame
C_body2gcrf = AttitudeTNRToGCRF( ...
        initial_satellite_position_gcrf, ...
        initial_satellite_velocity_gcrf);

facet_orientations_gcrf = C_body2gcrf * structure_parameters.facet_orientations_body;

% (2d) Compute the effective drag area
% The solar panel points in the direction of the sun.
solar_panel_facet_mask =strcmp(structure_parameters.facet_descriptions, "Solar");
facet_orientations_gcrf(:, solar_panel_facet_mask) = relative_sun_direction_gcrf;

facet_areas = structure_parameters.facet_areas;

ignore_facet_backwards = false(size(facet_areas));
ignore_facet_backwards(solar_panel_facet_mask) = true;

effective_area = ComputeEffectiveDragArea( ...
        initial_velocity_direction_gcrf, ...
        facet_orientations_gcrf, ...
        facet_areas, ...
        ignore_facet_backwards);

% (2e) Now compute the drag
drag_acceleration_gcrf = SimpleDragGCRF(...
    initial_satellite_position_gcrf, ...
    initial_satellite_velocity_gcrf, ...
    initial_coefficient_of_drag, ...
    effective_area, ...
    drag_parameters.mass, ...
    drag_parameters.reference_air_density, ...
    drag_parameters.reference_radial_distance, ...
    drag_parameters.decay_rate, ...
    drag_parameters.earth_rotation_rate);

% (3) Sun and Moon Gravity
moon_gravitational_acceleration_gcrf = ...
    SimpleTwoBodyOffsetGravitationalAccelerationGCRF(...
        (initial_satellite_position_gcrf - moon_position_gcrf), ...
        gravity_parameters.moon_gravitational_parameter);

sun_gravitational_acceleration_gcrf = ...
    SimpleTwoBodyOffsetGravitationalAccelerationGCRF(...
        (initial_satellite_position_gcrf - sun_position_gcrf), ...
        gravity_parameters.sun_gravitational_parameter);

% (4) Solar Radiation Pressure
warning("Solar radiation pressure not yet implemented!");

%% Measurements
% Transform the station position and velocity into the GCRF frame
% Vallado, algorithm 24
station_positions_gcrf  = C_itrf2gcrf * station_positions_itrf;
station_velocities_gcrf = C_mod2gcrf * C_tod2mod * C_pef2tod * (cross(repmat([0;0;omega], 1, 3), C_itrf2pef * station_positions_itrf, 1));

%% Predict the first measurement
measurements                     = load("provided/LEO_DATA_Apparent-1.mat");
measurements.station_number      = measurements.LEO_DATA_Apparent(:, 1);
measurements.time                = measurements.LEO_DATA_Apparent(:, 2);
measurements.apparent_range      = measurements.LEO_DATA_Apparent(:, 3) * Units.kilometers;
measurements.apparent_range_rate = measurements.LEO_DATA_Apparent(:, 4) * Units.kilometers / Units.seconds;
measurements.LEO_DATA_Apparent   = [];

satellite_position_at_rx_gcrf = initial_satellite_position_gcrf;
satellite_velocity_at_rx_gcrf = initial_satellite_velocity_gcrf;
station_position_at_rx_gcrf   = station_positions_gcrf(:, measurements.station_number(1));
station_velocity_at_rx_gcrf   = station_velocities_gcrf(:, measurements.station_number(1));

[apparent_range, apparent_range_rate, light_time] = ApparentRangeAndRangeRate( ...
        satellite_position_at_rx_gcrf, ...
        satellite_velocity_at_rx_gcrf, ...
        station_position_at_rx_gcrf, ...
        station_velocity_at_rx_gcrf, ...
        "threshold", 0.01 / Constants.SPEED_OF_LIGHT);

H = ApparentRangeAndRangeRateJacobian( ...
        satellite_position_at_rx_gcrf, ...
        satellite_velocity_at_rx_gcrf, ...
        station_position_at_rx_gcrf, ...
        station_velocity_at_rx_gcrf, ...
        light_time);
