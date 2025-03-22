function [y, Pyy, Pxy, H] = CubatureTransform(x, Pxx, h, R)
    % Applies the cubature transform to a state vector.
    %
    % Suppose a measurement is modeled as a nonlinear function of a state vector plus additive white Gaussin noise:
    %   y = h(x) + w
    % Here, x is the state, y is the measurement, h is the nonlinear function, and w is zero-mean Gaussian noise with
    % covariance R.
    %
    % The cubature transform computes the approximate mean measurement y, assosciated measurement-measurement and
    % state-measurement covariances Pyy and Pyx, and linearized measurement map H from a set of sample points.
    % 
    % The cubature transform could also be to evaluate a discrete-time station transition system of the form
    %   x[k+1] = f(x[k]) + v[k]
    % Here, x[k] and x[k+1] is the state at time t[k] and t[k+1], respectively, f is the nonlinear discrete-time state
    % transition function, and v is zero-mean Gaussian noise with covariance Q. 
    % 
    % Requires:
    % - x:
    % - Pxx:
    % - h:
    % - R:
    %
    % Returns:
    % - y:
    % - Pyy:
    % - Pxy:
    % - H:
    %
    % References:
    % - https://www.researchgate.net/profile/Haran-Arasaratnam/publication/278158963_Cubature_Kalman_Filtering_Theory_Applications/links/557cf6ad08ae26eada8ca08b/Cubature-Kalman-Filtering-Theory-Applications.pdf
    arguments(Input)
        x(:, 1) double {mustBeReal, mustBeFinite}
        Pxx(:, :) double {mustBeReal, mustBeFinite}
        h(1, 1) function_handle
        R(:, :) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        y(:, 1) double {mustBeReal, mustBeFinite}
        Pyy(:, :) double {mustBeReal, mustBeFinite}
        Pxy(:, :) double {mustBeReal, mustBeFinite}
        H(:, :) double {mustBeReal, mustBeFinite}
    end

    num_states = numel(x);
    num_measurements = size(R, 1);

    % Input argument validation
    assert(all(size(Pxx) == [num_states, num_states]));
    assert(all(size(R) == [num_measurements, num_measurements]));

    % Create sample points
    Pxx_lower         = chol(Pxx, 'lower');
    chi               = sqrt(num_states) * [Pxx_lower, -Pxx_lower] + x;
    chi_bar           = x;
    num_sample_points = size(chi, 2);

    % Apply the measurement function to the sample points.
    y = h(chi);

    % Compute the measurement statistics.
    y_bar = (1 / num_sample_points) * sum(y, 2);
    Pyy = (1 / num_sample_points) * (y - y_bar) * transpose(y - y_bar) + R;
    Pxy = (1 / num_sample_points) * (chi - chi_bar) * transpose(y - y_bar);
    
    % Compute the linearized measurement mapping
    H = transpose(Pxx \ Pxy);

    % Output
    y   = y_bar;
    Pyy = Pyy;
    Pxy = Pxy;
    H   = H;

    % Output argument validation
    assert(numel(y) == num_measurements);
    assert(all(size(Pyy) == [num_measurements, num_measurements]));
    assert(all(size(Pxy) == [num_states, num_measurements]));
    assert(all(size(H) == [num_measurements, num_states]));
end