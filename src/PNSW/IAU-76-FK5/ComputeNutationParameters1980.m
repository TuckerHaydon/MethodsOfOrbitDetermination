function [ ...
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
        obliquity_earth_orientation_parameter_deg)
    % Computes the nutation parameters using the IAU-76/FK5 1980 nutation model.
    % 
    % Requires:
    % - julian_centuries_since_epoch:
    %   - (1, 1) double.
    %   - The number of Julian centuries since epoch.
    % - moon_mean_anomaly_deg:
    %  - (:, :) double.
    %  - The mean anomaly of the moon.
    %  - Unit degrees.
    % - sun_mean_anomaly_deg:
    %   - (:, :) double.
    %   - The mean anomaly of the sun.
    %   - Unit degrees.
    % - moon_mean_argument_of_latitude_deg:
    %   - (:, :) double.
    %   - The mean argument of latitude of the moon.
    %   - Unit degrees.
    % - sun_mean_elongation_deg:
    %   - (:, :) double.
    %   - The mean elongation of the sun.
    %   - Unit degrees.
    % - moon_mean_longitude_of_ascending_node_deg:
    %   - (:, :) double.
    %   - The mean longitude of the ascending node of the moon.
    %   - Unit degrees.
    % - longitude_earth_orientation_parameter_deg:
    %   - (1, 1) double.
    %   - The earth orientation parameter adjustment for the longitudinal nutation.
    %   - Unit degrees.
    % - obliquity_earth_orientation_parameter_deg:
    %   - (1, 1) double.
    %   - The earth orientation parameter adjustment for the obliquity nutation.
    %   - Unit degrees.
    % 
    % Returns:
    % - nutation_in_longitude_1980_deg:
    %   - (1, 1) double.
    %   - The nutation in longitude.
    %   - Unit degrees.
    % - nutation_in_obliquity_1980_deg:
    %   - (1, 1) double.
    %   - The nutation is obliquity.
    %   - Unit degrees.
    % - true_obliquity_1980_deg:
    %  - (1, 1) double.
    %  - The true obliquity.
    %  - Unit degrees.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th edition
    arguments(Input)
        julian_centuries_since_epoch(1, 1) double {mustBeReal, mustBeFinite}
        moon_mean_anomaly_deg(1, 1) double {mustBeReal, mustBeFinite}
        sun_mean_anomaly_deg(1, 1) double {mustBeReal, mustBeFinite}
        moon_mean_argument_of_latitude_deg(1, 1) double {mustBeReal, mustBeFinite}
        sun_mean_elongation_deg(1, 1) double {mustBeReal, mustBeFinite}
        moon_mean_longitude_of_ascending_node_deg(1, 1) double {mustBeReal, mustBeFinite}
        mean_obliquity_of_ecliptic_1980_deg(1, 1) double {mustBeReal, mustBeFinite}
        longitude_earth_orientation_parameter_deg(1, 1) double {mustBeReal, mustBeFinite} = 0
        obliquity_earth_orientation_parameter_deg(1, 1) double {mustBeReal, mustBeFinite} = 0
    end

    arguments(Output)
        nutation_in_longitude_1980_deg(1, 1) double {mustBeReal, mustBeFinite}
        nutation_in_obliquity_1980_deg(1, 1) double {mustBeReal, mustBeFinite}
        true_obliquity_1980_deg(1, 1) double {mustBeReal, mustBeFinite}
    end

    % Rename
    TTT                          = julian_centuries_since_epoch;
    M_moon_deg                   = moon_mean_anomaly_deg;
    M_sun_deg                    = sun_mean_anomaly_deg;
    u_bar_moon_deg               = moon_mean_argument_of_latitude_deg;
    D_sun_deg                    = sun_mean_elongation_deg;
    lambda_bar_ecliptic_deg      = moon_mean_longitude_of_ascending_node_deg;
    epsilon_bar_1980_deg         = mean_obliquity_of_ecliptic_1980_deg;
    delta_delta_psi_1980_deg     = longitude_earth_orientation_parameter_deg;
    delta_delta_epsilon_1980_deg = obliquity_earth_orientation_parameter_deg;

    % Load the table of 106 nutation parameters
    nutation_data = LoadNutationData();
    a1 = nutation_data(:, 1);
    a2 = nutation_data(:, 2);
    a3 = nutation_data(:, 3);
    a4 = nutation_data(:, 4);
    a5 = nutation_data(:, 5);

    A = nutation_data(:, 6);
    B = nutation_data(:, 7);
    C = nutation_data(:, 8);
    D = nutation_data(:, 9);

    % What do we use this for?
    idx = nutation_data(:, 10);

    % Vallado 3-86
    ap_deg = a1 .* M_moon_deg + a2 .* M_sun_deg + a3 .* u_bar_moon_deg + a4 .* D_sun_deg + a5 .* lambda_bar_ecliptic_deg;

    % A, B, C, D are in 0.0001 arcseconds / century
    delta_psi_1980_arcseconds     = sum((A + B .* TTT) .* sind(ap_deg), 1) / 1e4;
    delta_epsilon_1980_arcseconds = sum((C + D .* TTT) .* cosd(ap_deg), 1) / 1e4;

    delta_psi_1980_deg     = delta_psi_1980_arcseconds     * Conversions.ARCSECONDS_TO_DEGREES;
    delta_epsilon_1980_deg = delta_epsilon_1980_arcseconds * Conversions.ARCSECONDS_TO_DEGREES;

    % Wrap
    delta_psi_1980_deg     = wrapTo180(delta_psi_1980_deg);
    delta_epsilon_1980_deg = wrapTo180(delta_epsilon_1980_deg);

    % Add EOP corrections
    delta_psi_1980_deg     = delta_psi_1980_deg     + delta_delta_psi_1980_deg;
    delta_epsilon_1980_deg = delta_epsilon_1980_deg + delta_delta_epsilon_1980_deg;

    % Wrap
    delta_psi_1980_deg     = wrapTo180(delta_psi_1980_deg);
    delta_epsilon_1980_deg = wrapTo180(delta_epsilon_1980_deg);

    % True obliquity of the ecliptic
    % Vallado 3-88
    epsilon_1980_deg = epsilon_bar_1980_deg + delta_epsilon_1980_deg;

    % Wrap to 360
    epsilon_1980_deg = wrapTo180(epsilon_1980_deg);
    
    % Output
    nutation_in_longitude_1980_deg = delta_psi_1980_deg;
    nutation_in_obliquity_1980_deg = delta_epsilon_1980_deg;
    true_obliquity_1980_deg        = epsilon_1980_deg;
end

function [nutation_data] = LoadNutationData()
    % First four columns are unitless.
    % Next four columns are unit 0.0001 arcseconds / Julian century.
    nutation_data = [
         0  0  0  0  1 -171996 -174.2 92025  8.9    1;
         0  0  2 -2  2  -13187   -1.6  5736 -3.1    9;
         0  0  2  0  2   -2274   -0.2   977 -0.5   31;
         0  0  0  0  2    2062    0.2  -895  0.5    2;
         0  1  0  0  0    1426   -3.4    54 -0.1   10;
         1  0  0  0  0     712    0.1    -7   0    32;
         0  1  2 -2  2    -517    1.2   224 -0.6   11;
         0  0  2  0  1    -386   -0.4   200   0    33;
         1  0  2  0  2    -301     0    129 -0.1   34;
         0 -1  2 -2  2     217   -0.5   -95  0.3   12;
         1  0  0 -2  0    -158     0     -1   0    35;
         0  0  2 -2  1     129    0.1   -70   0    13;
        -1  0  2  0  2     123     0    -53   0    36;
         1  0  0  0  1      63    0.1   -33   0    38;
         0  0  0  2  0      63     0     -2   0    37;
        -1  0  2  2  2     -59     0     26   0    40;
        -1  0  0  0  1     -58   -0.1    32   0    39;
         1  0  2  0  1     -51     0     27   0    41;
         2  0  0 -2  0      48     0      1   0    14;
        -2  0  2  0  1      46     0    -24   0     3;
         0  0  2  2  2     -38     0     16   0    42;
         2  0  2  0  2     -31     0     13   0    45;
         2  0  0  0  0      29     0     -1   0    43;
         1  0  2 -2  2      29     0    -12   0    44;
         0  0  2  0  0      26     0     -1   0    46;
         0  0  2 -2  0     -22     0      0   0    15;
        -1  0  2  0  1      21     0    -10   0    47;
         0  2  0  0  0      17   -0.1     0   0    16;
         0  2  2 -2  2     -16    0.1     7   0    18;
        -1  0  0  2  1      16     0     -8   0    48;
         0  1  0  0  1     -15     0      9   0    17;
         1  0  0 -2  1     -13     0      7   0    49;
         0 -1  0  0  1     -12     0      6   0    19;
         2  0 -2  0  0      11     0      0   0     4;
        -1  0  2  2  1     -10     0      5   0    50;
         1  0  2  2  2      -8     0      3   0    54;
         0 -1  2  0  2      -7     0      3   0    53;
         0  0  2  2  1      -7     0      3   0    58;
         1  1  0 -2  0      -7     0      0   0    51;
         0  1  2  0  2       7     0     -3   0    52;
        -2  0  0  2  1      -6     0      3   0    20;
         0  0  0  2  1      -6     0      3   0    57;
         2  0  2 -2  2       6     0     -3   0    56;
         1  0  0  2  0       6     0      0   0    55;
         1  0  2 -2  1       6     0     -3   0    58;
         0  0  0 -2  1      -5     0      3   0    60;
         0 -1  2 -2  1      -5     0      3   0    21;
         2  0  2  0  1      -5     0      3   0    62;
         1 -1  0  0  0       5     0      0   0    61;
         1  0  0 -1  0      -4     0      0   0    24;
         0  0  0  1  0      -4     0      0   0    65;
         0  1  0 -2  0      -4     0      0   0    63;
         1  0 -2  0  0       4     0      0   0    64;
         2  0  0 -2  1       4     0     -2   0    22;
         0  1  2 -2  1       4     0     -2   0    23;
         1  1  0  0  0      -3     0      0   0    66;
         1 -1  0 -1  0      -3     0      0   0     6;
        -1 -1  2  2  2      -3     0      1   0    69;
         0 -1  2  2  2      -3     0      1   0    72;
         1 -1  2  0  2      -3     0      1   0    68;
         3  0  2  0  2      -3     0      1   0    71;
        -2  0  2  0  2      -3     0      1   0     5;
         1  0  2  0  0       3     0      0   0    67;
        -1  0  2  4  2      -2     0      1   0    82;
         1  0  0  0  2      -2     0      1   0    76;
        -1  0  2 -2  1      -2     0      1   0    74;
         0 -2  2 -2  1      -2     0      1   0     7;
        -2  0  0  0  1      -2     0      1   0    70;
         2  0  0  0  1       2     0     -1   0    75;
         3  0  0  0  0       2     0      0   0    77;
         1  1  2  0  2       2     0     -1   0    73;
         0  0  2  1  2       2     0     -1   0    78;
         1  0  0  2  1      -1     0      0   0    91;
         1  0  2  2  1      -1     0      1   0    85;
         1  1  0 -2  1      -1     0      0   0   102;
         0  1  0  2  0      -1     0      0   0    99;
         0  1  2 -2  0      -1     0      0   0    30;
         0  1 -2  2  0      -1     0      0   0    27;
         1  0 -2  2  0      -1     0      0   0   103;
         1  0 -2 -2  0      -1     0      0   0   100;
         1  0  2 -2  0      -1     0      0   0    94;
         1  0  0 -4  0      -1     0      0   0    80;
         2  0  0 -4  0      -1     0      0   0    83;
         0  0  2  4  2      -1     0      0   0   105;
         0  0  2 -1  2      -1     0      0   0    98;
        -2  0  2  4  2      -1     0      1   0    86;
         2  0  2  2  2      -1     0      0   0    90;
         0 -1  2  0  1      -1     0      0   0   101;
         0  0 -2  0  1      -1     0      0   0    97;
         0  0  4 -2  2       1     0      0   0    92;
         0  1  0  0  2       1     0      0   0    28;
         1  1  2 -2  2       1     0     -1   0    84;
         3  0  2 -2  2       1     0      0   0    93;
        -2  0  2  2  2       1     0     -1   0    81;
        -1  0  0  0  2       1     0     -1   0    79;
         0  0 -2  2  1       1     0      0   0    26;
         0  1  2  0  1       1     0      0   0    95;
        -1  0  4  0  2       1     0      0   0    87;
         2  1  0 -2  0       1     0      0   0    25;
         2  0  0  2  0       1     0      0   0   104;
         2  0  2 -2  1       1     0     -1   0    89;
         2  0 -2  0  1       1     0      0   0     8;
         1 -1  0 -2  0       1     0      0   0    88;
        -1  0  0  1  1       1     0      0   0    29;
        -1 -1  0  2  1       1     0      0   0    96;
         0  1  0  1  0       1     0      0   0   106;  
    ];
end