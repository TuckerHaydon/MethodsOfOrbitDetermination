function [semimajor_axis] = SemimajorAxisFromSemiparameter( ...
        semiparameter, ...
        eccentricity)
    % Converts an orbital semiparameter to a semimajor axis. The semimajor axis is defined for all orbits except for
    % parabolic orbits. An error is thrown if a parabolic orbit is identified. 
    %
    % Requires:
    % - semiparameter:
    %   - (:, :) double.
    %   - The orbital semiparameter.
    %   - Unit meters.
    % - eccentricity:
    %   - (:, :) double.
    %   - The corresponding orbital eccentricity.
    %   - Unitless.
    % 
    % Returns:
    % - semimajor_axis:
    %   - (:, :) double.
    %   - The corresponding orbital semimajor axis.
    %   - Unit meters.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications
    arguments(Input)
        semiparameter(:, :) double {mustBeReal, mustBeFinite, mustBePositive}
        eccentricity(:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    arguments(Output)
        semimajor_axis(:, :) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    assert(all(size(semiparameter) == size(eccentricity)));

    % If the eccentricity is unity, the trajectory is a parabola and the semimajor axis is undefined. 
    if any(IsParabolicOrbit(eccentricity), 'all')
        error("SemimajorAxisFromSemiparameter: Parabolic orbit has no semimajor axis!");
    end

    % Vallado 1-19
    semimajor_axis = semiparameter ./ (1.0 - eccentricity.^2);

    % Output check
    assert(all(size(semimajor_axis) == size(semiparameter)));
end