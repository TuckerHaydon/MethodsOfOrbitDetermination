function [orbital_period] = OrbitalPeriodFromSemimajorAxis( ...
        semimajor_axis, ...
        eccentricity, ...
        gravitational_parameter)
    % Compute the period of an elliptical orbit from its semimajor axis.
    % 
    % Requires:
    % - semimajor_axis
    %   - (:, :) double.
    %   - Semimajor axis of the orbit.
    %   - Unit meters.
    % - eccentricity:
    %   - (:, :) double.
    %   - The orbit eccentricity.
    %   - Unitless.
    % - gravitational_parameter:
    %   - (1, 1) double.
    %   - Gravitational parameter of the primary.
    %   - Unit m^3 / s^2.
    % 
    % Returns:
    % - orbital_period:
    %   - (:, :) double.
    %   - The corresponding orbital period.
    % 
    % References:
    % - https://en.wikipedia.org/wiki/Orbital_period#Small_body_orbiting_a_central_body
    arguments(Input)
        semimajor_axis(:, :) double {mustBeReal, mustBeFinite, mustBePositive}
        eccentricity(:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        gravitational_parameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    arguments(Output)
        orbital_period(:, :) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    % Check that the orbit is circular or elliptical
    if any(~IsCircularOrEllipticalOrbit(eccentricity), 'all')
        error("OrbitalPeriodFromSemimajorAxis: Only circular or elliptical orbits have orbital periods!");
    end

    % Rename
    a   = semimajor_axis;
    mu  = gravitational_parameter;

    % From Wikipedia
    T = 2.0 .* pi .* sqrt(a.^3 ./ mu);

    % Output
    orbital_period = T;

    % Output check
    assert(all(size(orbital_period) == size(semimajor_axis)));
end