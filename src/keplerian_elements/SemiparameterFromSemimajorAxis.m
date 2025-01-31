function [semiparameter] = SemiparameterFromSemimajorAxis( ...
        semimajor_axis, ...
        eccentricity)
    % Converts an orbital semimajor axis to a semiparameter. The semimajor axis is defined for all orbits except for
    % parabolic orbits. An error is thrown if a parabolic orbit is identified. 
    %
    %
    % Requires:
    % - semimajor_axis:
    %   - (:, :) double.
    %   - The orbital semimajor axis.
    %   - Unit meters.
    % - eccentricity:
    %   - (:, :) double.
    %   - The corresponding orbital eccentricity.
    %   - Unitless.
    % 
    % Returns:
    % - semiparameter:
    %   - (:, :) double.
    %   - The corresponding orbital semiparameter.
    %   - Unit meters.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications
    arguments(Input)
        semimajor_axis(:, :) double {mustBeReal, mustBeFinite, mustBePositive}
        eccentricity(:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    arguments(Output)
        semiparameter(:, :) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    assert(all(size(semimajor_axis) == size(eccentricity)));

    % If the eccentricity is unity, the trajectory is a parabola and the semimajor axis is undefined. 
    if any(IsParabolicOrbit(eccentricity), 'all')
        error("SemiparameterFromSemimajorAxis: Parabolic orbit has no semimajor axis!");
    end

    % Vallado 1-19
    semiparameter = semimajor_axis .* (1.0 - eccentricity.^2);

    % Output check
    assert(all(size(semiparameter) == size(semimajor_axis)));
end