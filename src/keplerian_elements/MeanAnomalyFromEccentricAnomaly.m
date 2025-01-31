function [mean_anomaly] = MeanAnomalyFromEccentricAnomaly( ...
        eccentric_anomaly, ...
        eccentricity)
    % Computes the mean anomaly from the eccentric anomaly.
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
    % - mean_anomaly:
    %   - (:, :) double.
    %   - The orbit mean anomaly.
    %   - Unit radians.
    % 
    % References:
    % - https://en.wikipedia.org/wiki/Mean_anomaly#Mean_anomaly_at_epoch
    % - Vallado, Fundamentals of Astrodynamics and Applications
    arguments(Input)
        eccentric_anomaly(:, :) double {mustBeReal, mustBeFinite}
        eccentricity(:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    arguments(Output)
        mean_anomaly(:, :) double {mustBeReal, mustBeFinite}
    end

    assert(all(size(eccentric_anomaly) == size(eccentricity)));

    % Check that the orbit is circular or elliptical
    if any(~IsCircularOrEllipticalOrbit(eccentricity), 'all')
        error("MeanAnomalyFromEccentricAnomaly: Only circular or elliptical orbits have eccentric anomalies!");
    end

    % Rename
    E = eccentric_anomaly;
    e = eccentricity;

    % Vallado, 2-4
    M = E - e .* sin(E);

    % Clamp to [0, 2*pi]
    M(M < 0) = M(M < 0) + 2.0 .* pi;

    % Output
    mean_anomaly = M;

    % Output check
    assert(all(size(mean_anomaly) == size(eccentric_anomaly)));
end