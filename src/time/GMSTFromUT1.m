function [GMST_deg] = GMSTFromUT1(num_julian_centuries_of_ut1_time_since_j2000)
    % Computes the GMST from a provided UT1-referenced time.
    %
    % Note: Use the JulianCenturiesSince2000 function to compute the input.
    %
    % Requires:
    % - num_julian_centuries_of_ut1_time_since_j2000:
    %   - (:, :) double.
    %   - The number of Julian centuries of ut1 time since j2000.
    % 
    % Returns:
    % - GMST_deg:
    %   - (:, :) double.
    %   - The corresponding Greenwich mean sidereal time.
    %   - Unit degrees.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th edition
    arguments(Input)        
        num_julian_centuries_of_ut1_time_since_j2000(:, :) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        GMST_deg(:, :) double {mustBeReal, mustBeFinite}
    end

    % Rename
    T_UT1 = num_julian_centuries_of_ut1_time_since_j2000;

    % Vallado, 3-46
    GMST_sec = ...
        67310.54841 ...
        + (876600 * Conversions.HOURS_TO_SECONDS + 8640184.812866) .* T_UT1 ...
        + 0.093104 .* T_UT1.^2 ...
        - 6.2e-6 .* T_UT1.^3;

    GMST_deg = GMST_sec * Conversions.SECONDS_TO_DEGREES;

    % Wrap
    GMST_deg = wrapTo360(GMST_deg);

    % % Rename
    % T_UT1 = vpa(num_julian_centuries_of_ut1_time_since_j2000);
    % 
    % % Vallado, 3-46
    % GMST_sec = ...
    %     vpa(67310.54841) ...
    %     + (vpa(876600) * vpa(Conversions.HOURS_TO_SECONDS) + vpa(8640184.812866)) .* T_UT1 ...
    %     + vpa(0.093104) .* T_UT1.^2 ...
    %     - vpa(6.2e-6) .* T_UT1.^3;
    % 
    % GMST_deg = GMST_sec * vpa(Conversions.SECONDS_TO_DEGREES);
    % 
    % % Wrap
    % GMST_deg = wrapTo360(GMST_deg);
    % 
    % GMST_deg = double(GMST_deg);
end