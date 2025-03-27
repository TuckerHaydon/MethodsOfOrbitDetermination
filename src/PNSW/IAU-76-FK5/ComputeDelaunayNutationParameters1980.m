function [ ...
            moon_mean_anomaly_deg, ...
            sun_mean_anomaly_deg, ...
            moon_mean_argument_of_latitude_deg, ...
            sun_mean_elongation_deg, ...
            moon_mean_longitude_of_ascending_node_deg] = ...
        ComputeDelaunayNutationParameters1980(terrestrial_time_julian_centuries_since_epoch)
    % Computes the Delaunay nutation parameters using the IAU 1980 nutation model.
    % 
    % Requires:
    % - terrestrial_time_julian_centuries_since_epoch:
    %   - (:, :) double.
    %   - The number of terrestrial time julian centuries since the J2000 epoch.
    %
    % Returns:
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
    %
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th edition
    arguments(Input)
        terrestrial_time_julian_centuries_since_epoch(:, :) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        moon_mean_anomaly_deg(:, :) double {mustBeReal, mustBeFinite}
        sun_mean_anomaly_deg(:, :) double {mustBeReal, mustBeFinite}
        moon_mean_argument_of_latitude_deg(:, :) double {mustBeReal, mustBeFinite}
        sun_mean_elongation_deg(:, :) double {mustBeReal, mustBeFinite}
        moon_mean_longitude_of_ascending_node_deg(:, :) double {mustBeReal, mustBeFinite}
    end

    % Rename
    T_tt = terrestrial_time_julian_centuries_since_epoch;
    r   = 1 * Conversions.REVOLUTIONS_TO_DEGREES;
    
    % Vallado 3-85
    % Delaunay Parameters
    M_moon_deg               = 134.96298139 + (1325 .* r  + 198.8673981) .* T_tt + 0.0086972 .* T_tt.^2 + 1.78e-5 .* T_tt.^3;
    M_sun_deg                = 357.52772333 + (99   .* r  + 359.0503400) .* T_tt - 0.0001603 .* T_tt.^2 - 3.3e-6  .* T_tt.^3;
    u_bar_moon_deg           = 93.27191028  + (1342 .* r  + 82.0175381)  .* T_tt - 0.0036825 .* T_tt.^2 + 3.1e-6  .* T_tt.^3;
    D_sun_deg                = 297.85036306 + (1236 .* r  + 307.1114800) .* T_tt - 0.0019142 .* T_tt.^2 + 5.3e-6  .* T_tt.^3;
    lambda_bar_ecliptic_deg  = 125.04452222 - (5    .* r  + 134.1362608) .* T_tt + 0.0020708 .* T_tt.^2 + 2.2e-6  .* T_tt.^3;

    % Wrap to 360
    M_moon_deg              = wrapTo360(M_moon_deg);
    M_sun_deg               = wrapTo360(M_sun_deg);
    u_bar_moon_deg          = wrapTo360(u_bar_moon_deg);
    D_sun_deg               = wrapTo360(D_sun_deg);
    lambda_bar_ecliptic_deg = wrapTo360(lambda_bar_ecliptic_deg);

    % Rename
    moon_mean_anomaly_deg                         = M_moon_deg;
    sun_mean_anomaly_deg                          = M_sun_deg;
    moon_mean_argument_of_latitude_deg            = u_bar_moon_deg;
    sun_mean_elongation_deg                       = D_sun_deg;
    moon_mean_longitude_of_ascending_node_deg     = lambda_bar_ecliptic_deg;
end