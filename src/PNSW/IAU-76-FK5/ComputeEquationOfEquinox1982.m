function [equation_of_equinox_1982_deg] = ComputeEquationOfEquinox1982( ...
        nutation_in_longitude_1980_deg, ...
        mean_obliquity_of_ecliptic_1980_deg, ...
        moon_mean_longitude_of_ascending_node_deg, ...
        include_1994_resolution_correction)
    % Compute the Equation of the Equinoxes for with the IAU 1982 model.
    %
    % Requires:
    % - nutation_in_longitude_1980_deg:
    %   - (:, :) double.
    %   - The nutation in longitude computed with the 1980 nutation model.
    %   - Unit degrees.
    % - mean_obliquity_of_ecliptic_1980_deg:
    %   - (:, :) double.
    %   - The mean obliquity of the ecliptic computed with the 1980 nutation model.
    %   - Unit degrees.
    % - moon_mean_longitude_of_ascending_node_deg:
    %   - (:, :) double.
    %   - The mean longitude of the ascending node Delaunay parameter.
    %   - Unit degrees.
    % - include_1994_resolution_correction:
    %   - (1, 1) logical.
    %   - Whether to include the IAU 1994 resolution correction.
    %   - Defaults to true.
    %
    % Returns:
    % - equation_of_equinox_1982_deg:
    %   - (:, :) double.
    %   - The equation of equinoxes for 1982.
    %   - Unit degrees.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th edition
    arguments(Input)
        nutation_in_longitude_1980_deg(:, :) double {mustBeReal, mustBeFinite}
        mean_obliquity_of_ecliptic_1980_deg(:, :) double {mustBeReal, mustBeFinite}
        moon_mean_longitude_of_ascending_node_deg(:, :) double {mustBeReal, mustBeFinite}
        include_1994_resolution_correction(1, 1) logical = true
    end

    arguments(Output)
        equation_of_equinox_1982_deg(:, :) double {mustBeReal, mustBeFinite}
    end

    assert(all(size(nutation_in_longitude_1980_deg) == size(mean_obliquity_of_ecliptic_1980_deg)), ...
        "Input must be same size.");
    assert(all(size(nutation_in_longitude_1980_deg) == size(moon_mean_longitude_of_ascending_node_deg)), ...
        "Input must be same size.");

    % Rename
    delta_psi_1980_deg      = nutation_in_longitude_1980_deg;
    epsilon_bar_1980_deg    = mean_obliquity_of_ecliptic_1980_deg;
    lambda_bar_ecliptic_deg = moon_mean_longitude_of_ascending_node_deg;

    % Vallado, 3-82
    if include_1994_resolution_correction
        eq_eq_1982_deg = delta_psi_1980_deg .* cosd(epsilon_bar_1980_deg) + ...
            (0.00264  .* Conversions.ARCSECONDS_TO_DEGREES) .* sind(lambda_bar_ecliptic_deg) + ...
            (0.000063 .* Conversions.ARCSECONDS_TO_DEGREES) .* sind(2 .* lambda_bar_ecliptic_deg);
    else
        eq_eq_1982_deg = delta_psi_1980_deg .* cosd(epsilon_bar_1980_deg);
    end

    % Output
    equation_of_equinox_1982_deg = eq_eq_1982_deg;
end