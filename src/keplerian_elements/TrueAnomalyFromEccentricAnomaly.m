function [true_anomaly] = TrueAnomalyFromEccentricAnomaly( ...
        eccentric_anomaly, ...
        eccentricity)
    % Converts the eccentric anomaly of an elliptical orbit to the true anomaly. 
    % 
    % Note, cannot compute the eccentric anomaly for non-elliptical orbits.
    %
    % Requires:
    % - eccentric_anomaly:
    %   - (:, :) double.
    %   - The orbit eccentric anomaly.
    %   - Unit radians. 
    % - eccentricity:
    %   - (:, :) double.
    %   - The orbit eccentricity.
    %   - Unitless.
    % 
    % Returns:
    % - true_anomaly:
    %   - (:, :) double.
    %   - The orbit true anomaly.
    %   - Unit radians.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications
    arguments(Input)
        eccentric_anomaly(:, :) double {mustBeReal, mustBeFinite}
        eccentricity(:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    arguments(Output)
        true_anomaly(:, :) double {mustBeReal, mustBeFinite}
    end

    assert(all(size(eccentric_anomaly) == size(eccentricity)));

    % Check that the orbit is elliptical
    if any(~IsCircularOrEllipticalOrbit(eccentricity), 'all')
        error("TrueAnomalyFromEccentricAnomaly: Only circular or elliptical orbits have eccentric anomalies!");
    end

    % Rename
    E = eccentric_anomaly;
    e = eccentricity;

    % Vallado, Section 2.2.6, Algorithm 6
    nu = atan2(sqrt(1 - e.^2) .* sin(E), cos(E) - e);

    % Clamp to [0, 2*pi]
    nu(nu < 0) = nu(nu < 0) + 2.0 .* pi;

    % Output
    true_anomaly = nu;
end