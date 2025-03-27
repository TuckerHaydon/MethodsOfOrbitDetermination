function [sun_poisition_tod] = ComputeSunPositionTOD( ...
        ut1_time, ...
        tdb_time)
    % Computes the sun position vector in the mean of date (MOD) frame.
    %
    % Requires:
    % - ut1_time:
    %   - (1, 1) datetime.
    %   - The current ut1 datetime.
    % - tdb_time:
    %   - (1, 1) datetime.
    %   - The current tdb datetime.
    % 
    % Returns:
    % - sun_position_tod:
    %   - (3, 1) double.
    %   - The sun position vector in the TOD frame.
    %   - Unit meters.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th Edition
    arguments(Input)
        ut1_time(1, 1) datetime
        tdb_time(1, 1) datetime
    end
    
    arguments(Output)
        sun_poisition_tod(3, 1) double {mustBeReal, mustBeFinite}
    end

    % Timezone should be empty
    assert(isempty(ut1_time.TimeZone));
    assert(isempty(tdb_time.TimeZone));

    % Using Valldo Algorithm 29
    ut1_julian_date = JulianDate( ...
        ut1_time.Year, ...
        ut1_time.Month, ...
        ut1_time.Day, ...
        ut1_time.Hour, ...
        ut1_time.Minute, ...
        ut1_time.Second);

    tdb_julian_date = JulianDate( ...
        tdb_time.Year, ...
        tdb_time.Month, ...
        tdb_time.Day, ...
        tdb_time.Hour, ...
        tdb_time.Minute, ...
        tdb_time.Second);

    T_ut1 = JulianCenturiesSinceJ2000(ut1_julian_date);
    T_tdb = JulianCenturiesSinceJ2000(tdb_julian_date);

    lambda_bar_deg = 280.460 + 36000.771285 .* T_ut1;
    M_deg          = 357.528 + 35999.050957 .* T_tdb;

    lambda_bar_deg = wrapTo360(lambda_bar_deg);
    M_deg          = wrapTo360(M_deg);
    
    lambda_deg = lambda_bar_deg + 1.915 .* sind(M_deg) + 0.020 .* sind(2 .* M_deg);

    epsilon_deg = 23.439291 - 0.01461 .* T_tdb;

    r_mag_au = 1.00014 - 0.01671 .* cosd(M_deg) - 0.00014 .* cosd(2 .* M_deg);

    r_vec_au = [
        r_mag_au .* cosd(lambda_deg);
        r_mag_au .* cosd(epsilon_deg) .* sind(lambda_deg);
        r_mag_au .* sind(epsilon_deg) .* sind(lambda_deg);
    ];

    sun_poisition_tod = r_vec_au .* Constants.ASTRONOMICAL_UNIT;
end

% function [] = Meeus()
% 
    % - Meeus, Astronomical Algorithms

%     T_tdb = JulianCenturiesSinceJ2000(tdb_julian_date);
% 
%     L0_deg = 280.46645 + 36000.76983 .* T_tdb + 0.0003032 .* T_tdb.^2;
%     M_deg  = 357.52910 + 35999.05030 .* T_tdb - 0.0001559 .* T_tdb.^2 - 0.00000048 .* T_tdb.^3;
% 
%     L0_deg = wrapTo360(L0_deg);
%     M_deg  = wrapTo360(M_deg);
% 
%     e = 0.016708617 - 0.000042037 .* T_tdb - 0.0000001236 .* T_tdb.^2;
% 
%     C_deg = (1.914600 - 0.004817 .* T_tdb - 0.000014 .* T_tdb.^2) .* sind(M_deg) ...
%         + (0.019993 - 0.000101 .* T_tdb) .* sind(2 * M_deg) ...
%         + (0.000290) .* sind(3 * M_deg);
% 
%     O_dot_deg = L_deg + C_deg;
%     nu_deg    = M_deg + C_deg;
% 
%     % J2000 correction
%     O_dot_deg = O_dot_deg - 0.01397 .* (tdb_time.Year - 2000);
% 
%     R_au = 1.000001018 .* (1 - e.^2) ./ (1 + e .* cosd(nu_deg));
% 
%     Omega_deg = 125.04 - 1934.137 .* T_tdb;
%     lambda_deg = O_dot_deg - 0.00569 - 0.00478 .* sind(Omega_deg);
% 
%     epsilon_0_deg = 23 ...
%         + 26 .* Conversions.ARCMINUTES_TO_DEGREES ...
%         + 21.448 .* Conversions.ARCSECONDS_TO_DEGREES ...
%         - 46.8150 .* T_tdb     .* Conversions.ARCSECONDS_TO_DEGREES ...
%         + 0.00059 .* T_tdb.^2  .* Conversions.ARCSECONDS_TO_DEGREES ...
%         + 0.001813 .* T_tdb.^3 .* Conversions.ARCSECONDS_TO_DEGREES;
% 
%     epsilon_0_deg = epsilon_0_deg + ...
%         0.00256 .* cosd(Omega_deg);
% 
%     alpha_deg = atan2d(cosd(epsilon_0_deg) .* sind(lambda_deg), cosd(lambda_deg));
%     delta_deg = asind(sind(epsilon_0_deg) .* sind(lambda_deg));
% end