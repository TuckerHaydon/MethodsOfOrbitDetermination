function [moon_position_gcrf] = ComputeMoonPositionGCRF( ...
        tdb_time, ...
        earth_radius)
    % Computes the sun position vector in the mean of date (MOD) frame.
    %
    % Requires:
    % - tdb_time:
    %   - (1, 1) datetime.
    %   - The TDB datetime for the moon's position.
    % 
    % Returns:
    % - moon_position_gcrf:
    %   - (3, 1) double.
    %   - The moon position vector in the GCRF frame.
    %   - Unit meters.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th Edition
    arguments(Input)
        tdb_time(1, 1) datetime
        earth_radius(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end
    
    arguments(Output)
        moon_position_gcrf(3, 1) double {mustBeReal, mustBeFinite}
    end

    % Timezone should be empty
    assert(isempty(tdb_time.TimeZone));

    % Using Valldo Algorithm 31
    tdb_julian_date = JulianDate( ...
        tdb_time.Year, ...
        tdb_time.Month, ...
        tdb_time.Day, ...
        tdb_time.Hour, ...
        tdb_time.Minute, ...
        tdb_time.Second);

    T_tdb = JulianCenturiesSinceJ2000(tdb_julian_date);

    lambda_deg = 218.32 + 481267.8813 .* T_tdb ...
        + 6.29 .* sind(134.9 + 477198.85 .* T_tdb) ...
        - 1.27 .* sind(259.2 - 413335.38 .* T_tdb) ...
        + 0.66 .* sind(235.7 + 890534.23 .* T_tdb) ...
        + 0.21 .* sind(269.9 + 954397.70 .* T_tdb) ...
        - 0.19 .* sind(357.5 + 35999.05  .* T_tdb) ...
        - 0.11 .* sind(186.6 + 966404.05 .* T_tdb);
    lambda_deg = wrapTo360(lambda_deg);

    phi_deg = ...
        + 5.13 .* sind(93.3 + 483202.03 .* T_tdb) ...
        + 0.28 .* sind(228.2 + 960400.87 .* T_tdb) ...
        - 0.28 .* sind(318.3 + 6003.18   .* T_tdb) ...
        - 0.17 .* sind(217.6 - 407332.20 .* T_tdb);
    phi_deg = wrapTo180(phi_deg);

    P_deg = 0.9508 ...
        + 0.0518 .* cosd(134.9 + 477198.85 .* T_tdb) ...
        + 0.0095 .* cosd(259.2 - 413335.38 .* T_tdb) ...
        + 0.0078 .* cosd(235.7 + 890534.23 .* T_tdb) ...
        + 0.0028 .* cosd(269.9 + 954397.70 .* T_tdb);
    P_deg = wrapTo180(P_deg);

    epsilon_bar_deg = 23.439291 - 0.0130042 .* T_tdb - 1.64e-7 .* T_tdb.^2 + 5.04e-7 .* T_tdb.^3;
    epsilon_bar_deg = wrapTo180(epsilon_bar_deg);

    r = earth_radius ./ sind(P_deg);

    moon_position_gcrf = r .* [
        cosd(phi_deg) .* cosd(lambda_deg);
        cosd(epsilon_bar_deg) .* cosd(phi_deg) .* sind(lambda_deg) - sind(epsilon_bar_deg) .* sind(phi_deg);
        sind(epsilon_bar_deg) .* cosd(phi_deg) .* sind(lambda_deg) + cosd(epsilon_bar_deg) .* sind(phi_deg);
    ];
end