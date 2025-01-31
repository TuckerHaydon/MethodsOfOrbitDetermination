function [eccentric_anomaly] = EccentricAnomalyFromTrueAnomaly( ...
        true_anomaly, ...
        eccentricity)
    % Converts the true anomaly of an elliptical orbit to the eccentric anomaly. 
    % 
    % Note, cannot compute the eccentric anomaly for non-elliptical orbits.
    %
    % Requires:
    % - true_anomaly:
    %   - (:, :) double.
    %   - The orbit true anomaly.
    %   - Unit radians.
    % - eccentricity:
    %   - (:, :) double.
    %   - The orbit eccentricity.
    %   - Unitless.
    % 
    % Returns:
    % - eccentric_anomaly:
    %   - (:, :) double.
    %   - The orbit eccentric anomaly.
    %   - Unit radians. 
    % 
    % References:
    % - https://en.wikipedia.org/wiki/Eccentric_anomaly#From_the_true_anomaly
    % - Vallado, Fundamentals of Astrodynamics and Applications
    arguments(Input)
        true_anomaly(:, :) double {mustBeReal, mustBeFinite}
        eccentricity(:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    arguments(Output)
        eccentric_anomaly(:, :) double {mustBeReal, mustBeFinite}
    end

    assert(all(size(true_anomaly) == size(eccentricity)));

    % Check that the orbit is circular or elliptical
    if any(~IsCircularOrEllipticalOrbit(eccentricity), 'all')
        error("EccentricAnomalyFromTrueAnomaly: Only circular or elliptical orbits have eccentric anomalies!");
    end
    
    % Rename
    nu = true_anomaly;
    e  = eccentricity;

    % Vallado, 2-9
    E = atan2(sqrt(1 - e.^2) .* sin(nu), e + cos(nu));

    % Clamp to [0, 2*pi]
    E(E < 0) = E(E < 0) + 2.0 .* pi;

    % Output
    eccentric_anomaly = E;

    % Output check
    assert(all(size(eccentric_anomaly) == size(true_anomaly)));
end