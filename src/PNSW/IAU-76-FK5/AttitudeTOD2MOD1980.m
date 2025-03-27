function [C_tod2mod] = AttitudeTOD2MOD1980( ...
        mean_obliquity_of_ecliptic_1980_deg, ...
        nutation_in_longitude_1980_deg, ...
        true_obliquity_1980_deg)
    % Computes a coordinate transformation matrix that transforms a vector expressed in the the true of date (TOD) frame
    % into the same vector expressed in the mean of date (MOD) frame. This transform accounts for the nutation of the 
    % Earth's pole.
    % 
    % Requires:
    % - mean_obliquity_of_ecliptic_1980_deg:
    %   - (:, 1) double.
    %   - The mean obliquity of the ecliptic computed using the 1980 IAU nutation model.
    %   - Unit degrees.
    % - nutation_in_longitude_1980_deg:
    %   - (:, 1) double.
    %   - The nutation in longitude computed using the 1980 IAU nutation model.
    %   - Unit degrees.
    % - true_obliquity_1980_deg:
    %   - (:, 1) double.
    %   - The true obliquity of the ecliptic computed using the 1980 IAU nutation model.
    %   - Unit degrees.
    % 
    % Returns:
    % - C_tod2mod:
    %   - (3, 3, :) double.
    %   - The corresponding coordinate transform matrices from PEF to TOD.
    %
    % References:
    % - Vallado, Fundamental of Astrodynamics and applications, 5th edition
    arguments(Input)
        mean_obliquity_of_ecliptic_1980_deg(:, 1) double {mustBeReal, mustBeFinite}
        nutation_in_longitude_1980_deg(:, 1) double {mustBeReal, mustBeFinite}
        true_obliquity_1980_deg(:, 1) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        C_tod2mod(3, 3, :) double {mustBeReal, mustBeFinite}
    end

    assert(all(size(mean_obliquity_of_ecliptic_1980_deg) == size(nutation_in_longitude_1980_deg)), ...
        "Inputs must be same size.");
    assert(all(size(mean_obliquity_of_ecliptic_1980_deg) == size(true_obliquity_1980_deg)), ...
        "Inputs must be same size.");

    % Rename
    epsilon_bar_1980_deg = mean_obliquity_of_ecliptic_1980_deg;
    delta_psi_1980_deg   = nutation_in_longitude_1980_deg;
    epsilon_1980_deg     = true_obliquity_1980_deg;

    C_tod2mod = pagemtimes( ...
        pagemtimes(ValladoROT1(-epsilon_bar_1980_deg), ValladoROT3(delta_psi_1980_deg)), ...
        ValladoROT1(epsilon_1980_deg));
end