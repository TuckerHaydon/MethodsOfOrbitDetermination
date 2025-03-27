function [potential] = EGM96GravitationalPotential(...
        position_itrf, ...
        earth_gravitational_parameter, ...
        earth_radius, ...
        degree)
    % Computes the EGM spherical harmonic gravitational potential up to order 20.
    %
    % Requires:
    % - position_itrf:
    %   - (3, :) double.
    %   - Positions vectors at which to evaluate the spherical harmonic model.
    %   - Unit meters.
    %   - ITRF frame.
    % - earth_gravitational_parameter:
    %   - (1, 1) double.
    %   - The Earth gravitational parameters in unit m^3 / s^2. 
    %   - Must be 3986004.415E+8 or there is a model mismatch.
    % - earth_radius:
    %   - (1, 1) double.
    %   - The Earth radius in unit meters.
    %   - Must be 6378136.3, otherwise there is a model mismatch.
    % - degree:
    %   - (1, 1) double.
    %   - The degree of the spherical harmonic model to evaluate.
    %   - Max 20.
    % 
    % Returns:
    % - potential:
    %   - (:, 1) double.
    %   - The gravitational potential.
    %   - Unit m^2 / s^2.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th Edition
    % - https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/gravity-SphericalHarmonics.pdf
    % - https://cddis.nasa.gov/926/egm96/getit.html     
    arguments(Input)
        position_itrf(3, :) double {mustBeReal, mustBeFinite}
        earth_gravitational_parameter(1, 1) double {mustBeReal, mustBeFinite}
        earth_radius(1, 1) double {mustBeReal, mustBeFinite}
        degree(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    arguments(Output)
        potential(:, 1) double {mustBeReal, mustBeFinite}
    end

    % Load the table into memory
    persistent GM Re C S
    if isempty(GM) || isempty(Re) || isempty(C) || isempty(S)
        [GM, Re, C, S] = EGM96GravitationalModel20();
    end
    max_degree = size(C, 1) - 1;
    
    assert( ...
        degree <= max_degree, ...
        sprintf("Model only supports up to degree %d", max_degree));

    % Check that the provided coefficients are the same as the model coefficients.
    assert(abs(earth_gravitational_parameter - GM) < 1e-9, "Model mismatch!");
    assert(abs(earth_radius - Re) < 1e-9, "Model mismatch!");

    % Rename
    mu = GM;

    num_vectors = size(position_itrf, 2);

    r = vecnorm(position_itrf, 2, 1);
    r = reshape(r, [], 1);

    mu_over_r = mu ./ r;
    Re_over_r = Re ./ r;

    % Geocentric longitude
    lambda = atan2(position_itrf(2, :), position_itrf(1, :));
    lambda = reshape(lambda, [], 1);

    % sin(geocentric latitude)
    sin_phi = reshape(position_itrf(3, :), [], 1) ./ r;

    % Precompute the cos(m * lambda) and sin(m * lambda) terms.
    m = reshape(0:1:degree, 1, []);
    cos_m_lambda = cos(m .* lambda);
    sin_m_lambda = sin(m .* lambda);

    potential = zeros(num_vectors, 1);

    % Model computed with Vallado 8-19.
    for ell = 2:1:degree

        % The coefficients are normalized, so the normalization factor must be applied to the legendre polynomials as
        % well.
        m_ = reshape(0:1:ell, 1, []);
        normalization_factor = sqrt(...
            (2 - (m_ == 0)) .* ...
            (2 .* ell + 1) .* ...
            factorial(ell - m_) ./ ...
            factorial(ell + m_));

        % Note that the matlab legendre function is off by a factor of ((-1).^m_) from its definition in gravitational
        % spherical harmonic literature.
        P_bar_unnormalized = transpose(legendre(ell, sin_phi)) .* (-1).^m_;
        P_bar              = normalization_factor .* P_bar_unnormalized;

        row_potential = zeros(num_vectors, 1);

        for m = 0:1:ell
            row_potential = row_potential + ...
                P_bar(:, m + 1) .* (...
                    C(ell + 1, m + 1) .* cos_m_lambda(:, m + 1) + ...
                    S(ell + 1, m + 1) .* sin_m_lambda(:, m + 1));
        end

        row_potential = row_potential .* power(Re_over_r, ell);

        potential = potential + row_potential;
    end

    potential = mu_over_r .* (1 + potential);
end
