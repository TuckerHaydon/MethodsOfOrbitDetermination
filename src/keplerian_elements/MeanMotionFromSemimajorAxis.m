function [mean_motion] = MeanMotionFromSemimajorAxis( ...
        semimajor_axis, ...
        eccentricity, ...
        gravitational_parameter)
    % Computes the mean motion of an elliptical or circular orbit from its semimajor axis. 
    %
    % Requires:
    % - semimajor_axis:
    %   - (:, :) double.
    %   - The semimajor axis of the orbit.
    %   - Unit meters.
    % - eccentricity:
    %   - (:, :) double.
    %   - The orbit eccentricity.
    %   - Unitless.
    % - gravitational_parameter:
    %   - (1, 1) double.
    %   - The gravitation parameter of the primary body.
    %   - Units m^3 / s^2.
    % 
    % Returns:
    % - mean_motion:
    %   - (:, :) double.
    %   - The mean motion (mean angular speed) of the orbit.
    %   - Unit rad/s.
    % 
    % References:
    % - https://en.wikipedia.org/wiki/Orbital_period#Small_body_orbiting_a_central_body
    arguments(Input)
        semimajor_axis(:, :) double {mustBeReal, mustBeFinite, mustBePositive}
        eccentricity(:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        gravitational_parameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    arguments(Output)
        mean_motion(:, :) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    % Check that the orbit is circular or elliptical
    if any(~IsCircularOrEllipticalOrbit(eccentricity), 'all')
        error("MeanMotionFromSemimajorAxis: Only circular or elliptical orbits have a mean motion!");
    end

    orbital_period = OrbitalPeriodFromSemimajorAxis( ...
        semimajor_axis, ...
        eccentricity, ...
        gravitational_parameter);

    % Mean motion (angular speed) is 2 * pi / period.
    mean_motion = 2.0 .* pi ./ orbital_period;

    % Output check
    assert(all(size(mean_motion) == size(semimajor_axis)));
end