%% Preamble
clc; clear all; close all;

%% Provided
position_itrf = [
    -1033.4793830; 
     7901.2952754; 
     6380.3565958;
] .* Units.kilometers;

velocity_itrf = [
    -3.225636520; 
    -2.872451450; 
     5.531924446;
] .* Units.kilometers / Units.seconds;

% TODO: Compute TT_JD from provided time
julian_date                  = 2453101.828154745;
julian_centuries_since_epoch = TerrestrialTimeJulianCenturiesSinceJ2000(julian_date);
GMST_deg                     = 312.8098943;

% TODO: Compute these from online

% Do notuse these. The EOP corrections are *not* used for applications referencing the original IAU-76/FK5 theory.
longitude_earth_orientation_parameter_deg = -0.052195 * Conversions.ARCSECONDS_TO_DEGREES;
obliquity_earth_orientation_parameter_deg = -0.003875 * Conversions.ARCSECONDS_TO_DEGREES;

polar_motion_deg = [
    -0.140682;
     0.333309;
] * Conversions.ARCSECONDS_TO_DEGREES;

%% Compute Precession
% This is correct
[...
    zeta_deg, ...
    theta_deg, ...
    z_deg] = ...
ComputePrecessionParameters1976( ...
    julian_centuries_since_epoch);

C_mod2gcrf = AttitudeMOD2GCRF1976( ...
        zeta_deg, ...
        theta_deg, ...
        z_deg);

%% Compute Nutation














%% Compute wobble due to polar motion
% This is correct.
C_itrf2pef = AttitudeITRF2PEF1976(polar_motion_deg);

%% Compute the Delaunay nutation parameters
% These are correct.
[ ...
    moon_mean_anomaly_deg, ...
    sun_mean_anomaly_deg, ...
    moon_mean_argument_of_latitude_deg, ...
    sun_mean_elongation_deg, ...
    moon_mean_longitude_of_ascending_node_deg] = ...
ComputeDelaunayNutationParameters1980( ...
    julian_centuries_since_epoch);

% expected_moon_mean_anomaly_deg                     = 314.9118590;
% expected_sun_mean_anomaly_deg                      = 91.9379931;
% expected_moon_mean_argument_of_latitude_deg        = 169.0968272;
% expected_sun_mean_elongation_deg                   = 196.7518116;
% expected_moon_mean_longitude_of_ascending_node_deg = 42.6046140;
% 
% moon_mean_anomaly_relative_error = (moon_mean_anomaly_deg - expected_moon_mean_anomaly_deg) / expected_moon_mean_anomaly_deg;
% sun_mean_anomaly_relative_error = (sun_mean_anomaly_deg - expected_sun_mean_anomaly_deg) / expected_sun_mean_anomaly_deg;
% moon_mean_argument_of_latitude_relative_error = (moon_mean_argument_of_latitude_deg - expected_moon_mean_argument_of_latitude_deg) / expected_moon_mean_argument_of_latitude_deg;
% sun_mean_elongation_relative_error = (sun_mean_elongation_deg - expected_sun_mean_elongation_deg) / expected_sun_mean_elongation_deg;
% moon_mean_longitude_of_ascending_node_relative_error = (moon_mean_longitude_of_ascending_node_deg - expected_moon_mean_longitude_of_ascending_node_deg) / expected_moon_mean_longitude_of_ascending_node_deg;
% 
% fprintf("moon_mean_anomaly_relative_error relative error:                     %0.3e\n", moon_mean_anomaly_relative_error);
% fprintf("sun_mean_anomaly_relative_error relative error:                      %0.3e\n", sun_mean_anomaly_relative_error);
% fprintf("moon_mean_argument_of_latitude_relative_error relative error:        %0.3e\n", moon_mean_argument_of_latitude_relative_error);
% fprintf("sun_mean_elongation_relative_error relative error:                   %0.3e\n", sun_mean_elongation_relative_error);
% fprintf("moon_mean_longitude_of_ascending_node_relative_error relative error: %0.3e\n", moon_mean_longitude_of_ascending_node_relative_error);

%% Compute the mean obliquity of the ecliptic
% This is correct
[ ...
    mean_obliquity_of_ecliptic_1980_deg] = ...
ComputeMeanObliquityOfEcliptic1980( ...
    julian_centuries_since_epoch);

% expected_mean_obliquity_of_ecliptic_1980_deg = 23.4387368;
% mean_obliquity_of_ecliptic_relative_error = (mean_obliquity_of_ecliptic_1980_deg - expected_mean_obliquity_of_ecliptic_1980_deg) / expected_mean_obliquity_of_ecliptic_1980_deg;
% fprintf("Mean Obliquity Of Ecliptic Relative Error:                           %0.3e\n", ...
%     mean_obliquity_of_ecliptic_relative_error);

%% Compute nutation parameters
% I believe this is correct
[ ...
    nutation_in_longitude_1980_deg, ...
    nutation_in_obliquity_1980_deg, ...
    true_obliquity_1980_deg] = ...
ComputeNutationParameters1980(...
    julian_centuries_since_epoch, ...
    moon_mean_anomaly_deg, ...
    sun_mean_anomaly_deg, ...
    moon_mean_argument_of_latitude_deg, ...
    sun_mean_elongation_deg, ...
    moon_mean_longitude_of_ascending_node_deg, ...
    mean_obliquity_of_ecliptic_1980_deg, ...
    longitude_earth_orientation_parameter_deg, ...
    obliquity_earth_orientation_parameter_deg);


% expected_nutation_in_longitude_1980_deg = -0.0034108;
% expected_nutation_in_obliquity_1980_deg = 0.0020316;
% expected_true_obliquity_1980_deg        = 23.4407685;
% 
% nutation_in_longitude_1980_relative_error = (nutation_in_longitude_1980_deg - expected_nutation_in_longitude_1980_deg) / expected_nutation_in_longitude_1980_deg;
% nutation_in_obliquity_1980_relative_error = (nutation_in_obliquity_1980_deg - expected_nutation_in_obliquity_1980_deg) / expected_nutation_in_obliquity_1980_deg;
% true_obliquity_1980_relative_error = (true_obliquity_1980_deg - expected_true_obliquity_1980_deg) / expected_true_obliquity_1980_deg;
% 
% % Expect 5 digits of precision
% fprintf("nutation_in_longitude_1980_relative_error relative error:            %0.3e\n", nutation_in_longitude_1980_relative_error);
% fprintf("nutation_in_obliquity_1980_relative_error relative error:            %0.3e\n", nutation_in_obliquity_1980_relative_error);
% % Expect 9 digits of precision
% fprintf("true_obliquity_1980_relative_error relative error:                   %0.3e\n", true_obliquity_1980_relative_error);

%% Compute nutation rotation
% I think this is correct.
C_tod2mod  = AttitudeTOD2MOD1980( ...
        mean_obliquity_of_ecliptic_1980_deg, ...
        nutation_in_longitude_1980_deg, ...
        true_obliquity_1980_deg);

%% Compute Spin
% I think this is correct
include_1994_resolution_correction = true;

[...
    equation_of_equinox_1982_deg] = ...
ComputeEquationOfEquinox1982( ...
    nutation_in_longitude_1980_deg, ...
    mean_obliquity_of_ecliptic_1980_deg, ...
    moon_mean_longitude_of_ascending_node_deg, ...
    include_1994_resolution_correction);

% Off by the last two digits
GAST_deg = GASTFromGMST(GMST_deg, equation_of_equinox_1982_deg);
% GAST_deg = GMST_deg + nutation_in_longitude_1980_deg * cosd(true_obliquity_1980_deg);
% expected_GAST_deg = 312.8067654;

C_pef2tod  = AttitudePEF2TOD1976(GAST_deg);



% longitude_earth_orientation_parameter_deg
% obliquity_earth_orientation_parameter_deg

%% Transform position into GCRF frame
position_pef  = C_itrf2pef * position_itrf; % This is correct.
position_tod  = C_pef2tod  * position_pef;  % This is correct
position_mod  = C_tod2mod  * position_tod;  % This is off by one meter
position_gcrf = C_mod2gcrf * position_mod;

expected_position_pef = [
    -1033.4750313;
     7901.3055856;
     6380.3445327;
] * Units.kilometers;

expected_position_tod = [
    5094.5147804;
    6127.3664612;
    6380.3445328;
] * Units.kilometers;

expected_position_mod = [
    5094.0283745;
    6127.8708164;
    6380.2485164;
] * Units.kilometers;

expected_position_gcrf = [
    5102.5089529;
    6123.011399;
    6378.1369338;
] * Units.kilometers;

pef_error = (expected_position_pef - position_pef)
tod_error =(expected_position_tod - position_tod)
mod_error =(expected_position_mod - position_mod)
gcrf_error = (expected_position_gcrf - position_gcrf)

C_itrf2gcrf    = C_mod2gcrf * C_tod2mod * C_pef2tod * C_itrf2pef;
position_gcrf2 = C_itrf2gcrf * position_itrf;

%% Printout
% position_gcrf_km = position_gcrf / Units.kilometers