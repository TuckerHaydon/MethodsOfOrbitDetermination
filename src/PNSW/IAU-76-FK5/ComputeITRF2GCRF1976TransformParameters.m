function [C_itrf2gcrf, C_itrf2pef, C_pef2tod, C_tod2mod, C_mod2gcrf, R_dot, omega] = ...
    ComputeITRF2GCRF1976TransformParameters(...
        utc_time, ...
        ut1_utc_sec, ...
        polar_motion_deg, ...
        longitude_earth_orientation_parameter_deg, ...
        obliquity_earth_orientation_parameter_deg, ...
        length_of_day_earth_orientation_parameter_sec)
    % Computes the set of frame transformations to take a vector expressed in the International Terrestrial Reference
    % Frame (ITRF; ECEF) and express it in the Geocentric Celestial Reference Frame (GCRF; ECI). Uses the 1976 IAU-FK6
    % Theory of Nutation.
    %
    % Requires:
    % - utc_time:
    %   - (1, 1) datetime.
    %   - The UTC time assocated with the transform.
    % - ut1_utc_sec:
    %   - (1, 1) double.
    %   - The offset of UT1 from UTC.
    %   - Unit seconds.
    % - polar_motion_deg:
    %   - (2, 1) double.
    %   - The polar motion Earth orientation parameter.
    %   - Unit degrees.
    % - longitude_earth_orientation_parameter_deg:
    %   - (1, 1) double.
    %   - The nutation in longitude earth orientation parameter.
    %   - Unit degrees.
    % - obliquity_earth_orientation_parameter_deg:
    %   - (1, 1) double.
    %   - The nutation in obliquity earth orientation parameter.
    %   - Unit degrees.
    % - length_of_day_earth_orientation_parameter_sec:
    %   - (1, 1) double.
    %   - The length of day earth orientation parameter.
    %   - Unit seconds.
    % 
    % Returns:
    % - C_itrf2gcrf:
    %   - (3, 3) double.
    %   - The transform from the ITRF to the GCRF frame.
    % - C_itrf2pef:
    %   - (3, 3) double.
    %   - The transform from the ITRF to the PEF frame.
    % - C_pef2tod:
    %   - (3, 3) double.
    %   - The transform from the PEF to the TOD frame.
    % - C_tod2mod:
    %   - (3, 3) double.
    %   - The transform from the TOD to the MOD frame.
    % - C_mod2gcrf:
    %   - (3, 3) double.
    %   - The transform from the MOD to the GCRF frame.
    % - GAST_deg:
    %   - (3, 3) double.
    %   - The Greenwich Apparent Sidereal Time (GAST).
    %   - Unit degrees.
    % 
    % References: 
    % - Vallado, Fundamentals of astrodynamics and applications, 5th edition
    arguments(Input)
        utc_time(1, 1) datetime
        ut1_utc_sec(1, 1) double {mustBeReal, mustBeFinite}
        polar_motion_deg(2, 1) double {mustBeReal, mustBeFinite}
        longitude_earth_orientation_parameter_deg(1, 1) double {mustBeReal, mustBeFinite}
        obliquity_earth_orientation_parameter_deg(1, 1) double {mustBeReal, mustBeFinite}
        length_of_day_earth_orientation_parameter_sec(1, 1) double {mustBeReal, mustBeFinite}
    end

    % Convert time
    leapseconds_table     = leapseconds;
    leapseconds_table_idx = find(utc_time > leapseconds_table.Date, 1, 'last');
    leapseconds_at_epoch  = leapseconds_table.CumulativeAdjustment(leapseconds_table_idx);
    leapseconds_at_epoch  = ((2*(leapseconds_table.Type(leapseconds_table_idx) == "+")) - 1) * leapseconds_at_epoch;
    
    % Initial difference of 10 seconds plus leap seconds
    tai_utc_sec = 10 + seconds(leapseconds_at_epoch); 

    % Fixed by definition
    tt_tai_sec = 32.184; 
    
    % Compute different times
    ut1_time = utc_time + seconds(ut1_utc_sec);
    tai_time = utc_time + seconds(tai_utc_sec);
    tt_time  = tai_time + seconds(tt_tai_sec);
    
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
    
    % Compute Precession
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
    
    % Compute Nutation
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
    
    % Compute Spin
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
    
    % Compute Wobble
    % This is correct.
    C_itrf2pef = AttitudeITRF2PEF1976(polar_motion_deg);
    
    % Compose the frame transform
    C_itrf2gcrf = C_mod2gcrf * C_tod2mod * C_pef2tod * C_itrf2pef;

    % Compute the time derivative of this rotation
    % Vallado, 3-78
    LOD = length_of_day_earth_orientation_parameter_sec;
    omega = 7.292115146706979e-5 * (1 - LOD / 86400);

    R_dot = [
        -omega * sind(GAST_deg), -omega * cosd(GAST_deg), 0;
        +omega * cosd(GAST_deg), -omega * sind(GAST_deg), 0;
        0, 0, 0;
    ];
end