%% Preamble
clc; clear all; close all;

%% Intermediate values taken from vallado
% Use these to check
vallado.polar_motion_deg = [
    -0.140682;
     0.333309;
] * Conversions.ARCSECONDS_TO_DEGREES;

vallado.longitude_earth_orientation_parameter_deg = -0.052195 * Conversions.ARCSECONDS_TO_DEGREES;
vallado.obliquity_earth_orientation_parameter_deg = -0.003875 * Conversions.ARCSECONDS_TO_DEGREES;

vallado.tt_julian_date                              = 2453101.828154745;
vallado.GMST_deg                                    = 312.8098943;
vallado.num_julian_centuries_of_tt_time_since_j2000 = 0.0426236319;

%% Provided
% Provided ITRF position and velocity
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

% Provided UTC date
utc_year    = 2004;
utc_month   = 04;
utc_day     = 6;
utc_hour    = 7;
utc_minute  = 51;
utc_seconds = 28.386009;
utc_time    = datetime(utc_year, utc_month, utc_day, utc_hour, utc_minute, utc_seconds);

%% Determine the Earth orientation parameters of the provided date
% Downloaded the parameters from https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
% FILE: finals.all.iau1980.txt
% There are 185 characters in a full line.
eop_file = "finals.all.iau1980.txt";

% Read all lines
eop_lines = readlines(eop_file);

% Strip trailing whitespace
eop_lines = strip(eop_lines, 'right');

% Remove any incomplete lines
num_expected_characters = 185;
eop_lines = eop_lines(strlength(eop_lines) == 185);

% Turn it into a character array
eop_lines = char(eop_lines);

%% Parse the EOP
% I believe error is the standard deviation.
% Col.#    Format  Quantity
% -------  ------  -------------------------------------------------------------
% 1-2      I2      year (to get true calendar year, add 1900 for MJD<=51543 or add 2000 for MJD>=51544)
% 3-4      I2      month number
% 5-6      I2      day of month
% 7        X       [blank]
% 8-15     F8.2    fractional Modified Julian Date (MJD UTC)
% 16       X       [blank]
% 17       A1      IERS (I) or Prediction (P) flag for Bull. A polar motion values
% 18       X       [blank]
% 19-27    F9.6    Bull. A PM-x (sec. of arc)
% 28-36    F9.6    error in PM-x (sec. of arc)
% 37       X       [blank]
% 38-46    F9.6    Bull. A PM-y (sec. of arc)
% 47-55    F9.6    error in PM-y (sec. of arc)
% 56-57    2X      [blanks]
% 58       A1      IERS (I) or Prediction (P) flag for Bull. A UT1-UTC values
% 59-68    F10.7   Bull. A UT1-UTC (sec. of time)
% 69-78    F10.7   error in UT1-UTC (sec. of time)
% 79       X       [blank]
% 80-86    F7.4    Bull. A LOD (msec. of time) -- NOT ALWAYS FILLED
% 87-93    F7.4    error in LOD (msec. of time) -- NOT ALWAYS FILLED
% 94-95    2X      [blanks]
% 96       A1      IERS (I) or Prediction (P) flag for Bull. A nutation values
% 97       X       [blank]
% 98-106   F9.3    Bull. A dPSI (msec. of arc)
% 107-115  F9.3    error in dPSI (msec. of arc)
% 116      X       [blank]
% 117-125  F9.3    Bull. A dEPSILON (msec. of arc)
% 126-134  F9.3    error in dEPSILON (msec. of arc)
% 135-144  F10.6   Bull. B PM-x (sec. of arc)
% 145-154  F10.6   Bull. B PM-y (sec. of arc)
% 155-165  F11.7   Bull. B UT1-UTC (sec. of time)
% 166-175  F10.3   Bull. B dPSI (msec. of arc)
% 176-185  F10.3   Bull. B dEPSILON (msec. of arc)
eop.daily.year                                        = int32(str2num(eop_lines(:, 1:2)));
eop.daily.month                                       = int32(str2num(eop_lines(:, 3:4)));
eop.daily.day                                         = int32(str2num(eop_lines(:, 5:6)));
% line 7 blank
eop.daily.mjd_utc                                     = str2double(string(eop_lines(:, 8:15)));
% line 16 blank
eop.daily.polar_motion_type                           = eop_lines(:, 17);
% line 18 blank
eop.daily.polar_motion_x_arcseconds                   = str2double(string(eop_lines(:, 19:27)));
eop.daily.polar_motion_x_arcseconds_error             = str2double(string(eop_lines(:, 28:36)));
% line 37 blank
eop.daily.polar_motion_y_arcseconds                   = str2double(string(eop_lines(:, 38:46)));
eop.daily.polar_motion_y_arcseconds_error             = str2double(string(eop_lines(:, 47:55)));
% lines 56, 57 blank
eop.daily.ut1_utc_type                                = eop_lines(:, 58);
eop.daily.ut1_utc_sec                                 = str2double(string(eop_lines(:, 59:68)));
eop.daily.ut1_utc_sec_error                           = str2double(string(eop_lines(:, 69:78)));
% line 79 blank
eop.daily.lod_millisec                                = str2double(string(eop_lines(:, 80:86)));
eop.daily.lod_millisec_error                          = str2double(string(eop_lines(:, 87:93)));
% lines 94, 95 blank
eop.daily.nutation_type                               = eop_lines(:, 96);
% line 97 blank
eop.daily.nutation_dPsi_milliarcsecond                = str2double(string(eop_lines(:, 98:106)));
eop.daily.nutation_dPsi_milliarcsecond_error          = str2double(string(eop_lines(:, 107:115)));
% line 116 blank
eop.daily.nutation_dEpsilon_milliarcsecond            = str2double(string(eop_lines(:, 117:125)));
eop.daily.nutation_dEpsilon_milliarcsecond_error      = str2double(string(eop_lines(:, 126:134)));
eop.monthly.polar_motion_x_arcseconds                 = str2double(string(eop_lines(:, 135:144)));
eop.monthly.polar_motion_y_arcseconds                 = str2double(string(eop_lines(:, 145:154)));
eop.monthly.ut1_utc_sec                               = str2double(string(eop_lines(:, 155:165)));
eop.monthly.nutation_dPsi_milliarcsecond              = str2double(string(eop_lines(:, 166:175)));
eop.monthly.nutation_dEpsilon_milliarcsecond          = str2double(string(eop_lines(:, 176:185)));

% Convert year into four digits
mask_1900 = eop.daily.mjd_utc <= 51543;
mask_2000 = ~mask_1900;
eop.daily.year(mask_1900) = eop.daily.year(mask_1900) + 1900;
eop.daily.year(mask_2000) = eop.daily.year(mask_2000) + 2000;

%% Find the data of interest
year_mask  = eop.daily.year  == utc_year;
month_mask = eop.daily.month == utc_month;
day_mask   = eop.daily.day   == utc_day;

data_mask  = year_mask & month_mask & day_mask;

assert(sum(data_mask) == 1, "More or less than 1 record found!");

data_idx = find(data_mask, 1, 'first');

polar_motion_arcseconds = [
    eop.daily.polar_motion_x_arcseconds(data_idx); 
    eop.daily.polar_motion_y_arcseconds(data_idx); 
];

ut1_utc_sec           = eop.daily.ut1_utc_sec(data_idx);
lod_ms                = eop.daily.lod_millisec(data_idx);
dPsi_milli_arcseconds = eop.daily.nutation_dPsi_milliarcsecond(data_idx);
dEps_milli_arcseconds = eop.daily.nutation_dEpsilon_milliarcsecond(data_idx);

polar_motion_deg                          = polar_motion_arcseconds * Conversions.ARCSECONDS_TO_DEGREES;
lod_sec                                   = lod_ms / 1000;
longitude_earth_orientation_parameter_deg = dPsi_milli_arcseconds / 1000 * Conversions.ARCSECONDS_TO_DEGREES;
obliquity_earth_orientation_parameter_deg = dEps_milli_arcseconds / 1000 * Conversions.ARCSECONDS_TO_DEGREES;

%% Convert time
leapseconds_table     = leapseconds;
leapseconds_table_idx = find(utc_time > leapseconds_table.Date, 1, 'last');
leapseconds_at_epoch  = leapseconds_table.CumulativeAdjustment(leapseconds_table_idx);
leapseconds_at_epoch  = ((2*(leapseconds_table.Type(leapseconds_table_idx) == "+")) - 1) * leapseconds_at_epoch;
tai_utc_sec           = 10 + seconds(leapseconds_at_epoch); % initial difference of 10 seconds plus leap seconds
tt_tai_sec            = 32.184; % fixed by definition
gps_tai_sec           = -19; % fixed by definition

% Compute Times
ut1_time = utc_time + seconds(ut1_utc_sec);
tai_time = utc_time + seconds(tai_utc_sec);
tt_time  = tai_time + seconds(tt_tai_sec);
% gps_time = tai_time + seconds(gps_tai_sec);

[ut1_julian_date, ~] = JulianDate( ...
    year(ut1_time), ...
    month(ut1_time), ...
    day(ut1_time), ...
    hour(ut1_time), ...
    minute(ut1_time), ...
    second(ut1_time));

[tt_julian_date, ~] = JulianDate( ...
    year(tt_time), ...
    month(tt_time), ...
    day(tt_time), ...
    hour(tt_time), ...
    minute(tt_time), ...
    second(tt_time));

% Compute the GMST. This requires the number of UT1 "centuries" since the J2000.0 epoch.
num_julian_centuries_of_ut1_time_since_j2000 = JulianCenturiesSinceJ2000(ut1_julian_date);
GMST_deg = GMSTFromUT1(num_julian_centuries_of_ut1_time_since_j2000);

% Compute the number of terrestrial time julian "centuries" since epoch. This is required for most of the following
% formulas.
num_julian_centuries_of_tt_time_since_j2000 = JulianCenturiesSinceJ2000(tt_julian_date);

%% Compute Precession
% This is correct.
[...
    zeta_deg, ...
    theta_deg, ...
    z_deg] = ...
ComputePrecessionParameters1976( ...
    num_julian_centuries_of_tt_time_since_j2000);

% Compute rotation
% This is correct.
C_mod2gcrf = AttitudeMOD2GCRF1976( ...
        zeta_deg, ...
        theta_deg, ...
        z_deg);

%% Compute Nutation
% Compute the Delaunay nutation parameters
% These are correct.
[ ...
    moon_mean_anomaly_deg, ...
    sun_mean_anomaly_deg, ...
    moon_mean_argument_of_latitude_deg, ...
    sun_mean_elongation_deg, ...
    moon_mean_longitude_of_ascending_node_deg] = ...
ComputeDelaunayNutationParameters1980( ...
    num_julian_centuries_of_tt_time_since_j2000);

% Compute the mean obliquity of the ecliptic
% This is correct.
[ ...
    mean_obliquity_of_ecliptic_1980_deg] = ...
ComputeMeanObliquityOfEcliptic1980( ...
    num_julian_centuries_of_tt_time_since_j2000);

% Compute nutation parameters
% This is correct.
[ ...
    nutation_in_longitude_1980_deg, ...
    nutation_in_obliquity_1980_deg, ...
    true_obliquity_1980_deg] = ...
ComputeNutationParameters1980(...
    num_julian_centuries_of_tt_time_since_j2000, ...
    moon_mean_anomaly_deg, ...
    sun_mean_anomaly_deg, ...
    moon_mean_argument_of_latitude_deg, ...
    sun_mean_elongation_deg, ...
    moon_mean_longitude_of_ascending_node_deg, ...
    mean_obliquity_of_ecliptic_1980_deg, ...
    longitude_earth_orientation_parameter_deg, ...
    obliquity_earth_orientation_parameter_deg);

% Compute nutation rotation
% This is correct.
C_tod2mod  = AttitudeTOD2MOD1980( ...
        mean_obliquity_of_ecliptic_1980_deg, ...
        nutation_in_longitude_1980_deg, ...
        true_obliquity_1980_deg);

%% Compute Spin
% Correction is applied for dates after Feb 27, 1997.
include_1994_resolution_correction = utc_time > datetime(1997, 2, 27, 0, 0, 0);

% Compute the equation of equinoxes.
% This is correct, I think...
[...
    equation_of_equinox_1982_deg] = ...
ComputeEquationOfEquinox1982( ...
    nutation_in_longitude_1980_deg, ...
    mean_obliquity_of_ecliptic_1980_deg, ...
    moon_mean_longitude_of_ascending_node_deg, ...
    include_1994_resolution_correction);

% Compute the GAST.
% This is correct.
GAST_deg = GASTFromGMST(GMST_deg, equation_of_equinox_1982_deg);

% Compute the rotation.
% This is correct.
C_pef2tod  = AttitudePEF2TOD1976(GAST_deg);

%% Compute Wobble
% This is correct.
C_itrf2pef = AttitudeITRF2PEF1976(polar_motion_deg);

%% Earth precession term
LOD = lod_ms / 1000;
omega = 7.292115146706979e-5 * (1 - LOD / 86400);

R_dot = [
    -omega * sind(GAST_deg), -omega * cosd(GAST_deg), 0;
    +omega * cosd(GAST_deg), -omega * sind(GAST_deg), 0;
    0, 0, 0;
];

%% Transform position into GCRF frame
position_pef  = C_itrf2pef * position_itrf;
position_tod  = C_pef2tod  * position_pef;
position_mod  = C_tod2mod  * position_tod;
position_gcrf = C_mod2gcrf * position_mod;

velocity_pef  = C_itrf2pef * velocity_itrf;
velocity_tod  = C_pef2tod  * velocity_pef + R_dot * position_pef;
velocity_mod  = C_tod2mod  * velocity_tod;
velocity_gcrf = C_mod2gcrf * velocity_mod;

%% Compare against book answers
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

expected_velocity_gcrf = [
    -4.74322016;
    0.79053650;
    5.533755272;
] * Units.kilometers / Units.seconds;

% TOD is wrong, because its computation is inconsistent in the book.
pef_error  = (expected_position_pef  - position_pef)
tod_error  = (expected_position_tod  - position_tod)
mod_error  = (expected_position_mod  - position_mod)
gcrf_error = (expected_position_gcrf - position_gcrf)

gcrf_velocity_error = (expected_velocity_gcrf - velocity_gcrf)


diff_polar_motion = [polar_motion_deg, vallado.polar_motion_deg, polar_motion_deg - vallado.polar_motion_deg]
diff_eop = [
    longitude_earth_orientation_parameter_deg, vallado.longitude_earth_orientation_parameter_deg, longitude_earth_orientation_parameter_deg - vallado.longitude_earth_orientation_parameter_deg;
    obliquity_earth_orientation_parameter_deg, vallado.obliquity_earth_orientation_parameter_deg, obliquity_earth_orientation_parameter_deg - vallado.obliquity_earth_orientation_parameter_deg;
]

diff_tt_julian_date = [
    tt_julian_date, vallado.tt_julian_date, tt_julian_date - vallado.tt_julian_date
]

diff_GMST_deg = [
    GMST_deg, vallado.GMST_deg, GMST_deg - vallado.GMST_deg;
]

diff_num_julian_centuries_of_tt_time_since_j2000 = [
    num_julian_centuries_of_tt_time_since_j2000, vallado.num_julian_centuries_of_tt_time_since_j2000, num_julian_centuries_of_tt_time_since_j2000 - vallado.num_julian_centuries_of_tt_time_since_j2000
]
