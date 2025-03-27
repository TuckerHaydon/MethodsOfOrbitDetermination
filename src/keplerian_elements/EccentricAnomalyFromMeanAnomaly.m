function [eccentric_anomaly] = EccentricAnomalyFromMeanAnomaly( ...
        mean_anomaly, ...
        eccentricity, ...
        options)
    % Computes the eccentric anomaly from the mean anomaly.
    %
    % Note, this function is not vectorized due to an iterative loop.
    % 
    % Requires:
    % - mean_anomaly:
    %   - (1, 1) double.
    %   - The orbit mean anomaly.
    %   - Unit radians.
    % - eccentricity:
    %   - (1, 1) double.
    %   - The orbit eccentricity.
    %   - Unitless.
    % - options (struct):
    %   - abs_tolerance:
    %     - (1, 1) double. 
    %     - The absolute tolerance at which to stop the iteration.
    %   - rel_tolerance:
    %     - (1, 1) double. 
    %     - The relative tolerance at which to stop the iteration.
    %   - max_iterations:
    %     - (1, 1) double. 
    %     - The maximum number of iterations.
    %
    % Returns:
    % - eccentric_anomaly:
    %   - (1, 1) double.
    %   - The orbit eccentric anomaly.
    %   - Unit radians.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, Algorithm 2
    arguments(Input)
        mean_anomaly(1, 1) double {mustBeReal, mustBeFinite, mustBeIn2Pi}
        eccentricity(1, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        options.abs_tolerance(1, 1) double {mustBeReal, mustBeFinite, mustBePositive} = 1e-12
        options.rel_tolerance(1, 1) double {mustBeReal, mustBeFinite, mustBePositive} = 1e-12
        options.max_iterations(1, 1) double {mustBeReal, mustBeFinite, mustBePositive} = 30
    end

    arguments(Output)
        eccentric_anomaly(1, 1) double {mustBeReal, mustBeFinite}
    end

    % Check that the orbit is circular or elliptical
    if any(~IsCircularOrEllipticalOrbit(eccentricity), 'all')
        error("EccentricAnomalyFromMeanAnomaly: Only circular or elliptical orbits have eccentric anomalies!");
    end
    
    % Rename
    M = mean_anomaly;
    e = eccentricity;

    if (-pi < M && M < 0) || (pi < M && M < 2*pi)
        E = M - e;
    else
        E = M + e;
    end

    % Newton-Raphson iteration
    abs_tolerance  = options.abs_tolerance;
    rel_tolerance  = options.rel_tolerance;
    max_iterations = options.max_iterations;

    iteration = 0;
    E_k       = E;

    while true
        E_kp1 = E_k + (M - E_k + e .* sin(E_k)) ./ (1.0 - e .* cos(E_k));

        iteration = iteration + 1;

        % Tolerance break condition
        abs_error = abs(E_kp1 - E_k);
        rel_error = abs_error / E_k;
        if abs_error < abs_tolerance || rel_error < rel_tolerance
            E = E_kp1;
            break;
        end

        % Iteration break condition
        if iteration == max_iterations
            warn("EccentricAnomalyFromMeanAnomaly: Max iterations exceeded (%d)!", iteration);
            E = E_kp1;
            break;
        end

        % Else continue
        E_k = E_kp1;
    end

    % Output 
    eccentric_anomaly = E;
end

function [bool] = mustBeIn2Pi(x)
    bool = abs(x) <= 2 * pi;
end