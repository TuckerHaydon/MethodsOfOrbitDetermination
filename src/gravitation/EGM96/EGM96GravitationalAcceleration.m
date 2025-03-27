function [acceleration_itrf] = EGM96GravitationalAcceleration(...
        position_itrf, ...
        earth_gravitational_parameter, ...
        earth_radius, ...
        degree)
    % Computes the EGM spherical harmonic gravitational acceleration with a model order up to 20.
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
    % - acceleration_itrf:
    %   - (3, :) double.
    %   - The gravitational potential acceleration.
    %   - Unit m / s^2.
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
        acceleration_itrf(3, :) double {mustBeReal, mustBeFinite}
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

    % Unnormalize the table
    persistent C_unnormalized S_unnormalized
    if isempty(C_unnormalized) || isempty(S_unnormalized)
        [C_unnormalized, S_unnormalized] = UnnormalizeCoefficients(C, S);
    end

    % Rename
    mu = GM;

    num_vectors = size(position_itrf, 2);

    r = vecnorm(position_itrf, 2, 1);
    r = reshape(r, [], 1);

    mu_over_r2 = mu ./ (r.^2);
    Re_over_r  = Re ./ r;

    % Geocentric longitude
    lambda = atan2(position_itrf(2, :), position_itrf(1, :));
    lambda = reshape(lambda, [], 1);

    % Geocentric latitude
    phi = asin(reshape(position_itrf(3, :), [], 1) ./ r);

    sin_phi = sin(phi);
    cos_phi = cos(phi);

    cos_lambda = cos(lambda);
    sin_lambda = sin(lambda);

    % We need to compute dP(sin(phi)) / dphi later.
    dphi = 1e-3;
    sin_phi_plus  = sin(phi + dphi);
    sin_phi_minus = sin(phi - dphi);

    % Precompute the cos(m * lambda) and sin(m * lambda) terms.
    m = reshape(0:1:degree, 1, []);
    cos_m_lambda = cos(m .* lambda);
    sin_m_lambda = sin(m .* lambda);
    clear m

    a_r      = zeros(num_vectors, 1);
    a_phi    = zeros(num_vectors, 1);
    a_lambda = zeros(num_vectors, 1);

    % Following the method shown in
    % https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/gravity-SphericalHarmonics.pdf
    for ell = 0:1:degree

        % Note that the matlab legendre function is off by a factor of ((-1).^m_) from its definition in gravitational
        % spherical harmonic literature.
        m_ = reshape(0:1:ell, 1, []);
        P = transpose(legendre(ell, sin_phi)) .* (-1).^m_;

        % Computing dP(sin(phi)) / dphi with finite difference.
        P_plus  = transpose(legendre(ell, sin_phi_plus))  .* (-1).^m_;
        P_minus = transpose(legendre(ell, sin_phi_minus)) .* (-1).^m_;
        dP_dphi = (P_plus - P_minus) ./ (2 * dphi);

        a_r_row      = zeros(num_vectors, 1);
        a_phi_row    = zeros(num_vectors, 1);
        a_lambda_row = zeros(num_vectors, 1);

        for m = 0:1:ell

            C_ell_m = C_unnormalized(ell + 1, m + 1);
            S_ell_m = S_unnormalized(ell + 1, m + 1);
            P_ell_m = P(:, m + 1);
            
            a_r_row = a_r_row + ...
                 P_ell_m .* (C_ell_m .* cos_m_lambda(:, m + 1) + S_ell_m .* sin_m_lambda(:, m + 1));

            % These two loops start at 1, not zero.
            if 0 ~= ell
                a_phi_row = a_phi_row + ...
                    dP_dphi(:, m + 1) .* (C_ell_m .* cos_m_lambda(:, m + 1) + S_ell_m .* sin_m_lambda(:, m + 1));
    
                a_lambda_row = a_lambda_row + ...
                    m .* P_ell_m ./ cos_phi .* (S_ell_m .* cos_m_lambda(:, m + 1) - C_ell_m .* sin_m_lambda(:, m + 1));
            end
        end

        Re_over_r_to_ell = power(Re_over_r, ell);

        a_r_row = Re_over_r_to_ell .* (ell + 1) .* a_r_row;
        a_r = a_r + a_r_row;

        a_phi_row = Re_over_r_to_ell .* a_phi_row;
        a_phi = a_phi + a_phi_row;

        a_lambda_row = Re_over_r_to_ell .* a_lambda_row;
        a_lambda = a_lambda + a_lambda_row;
    end

    a_r      = -mu_over_r2 .* a_r;
    a_phi    = +mu_over_r2 .* a_phi;
    a_lambda = +mu_over_r2 .* a_lambda;

    % This is the transform matrix.
    % T = [
    %     cos_phi .* cos_lambda,    -sin_phi .* cos_lambda,    -sin_lambda;
    %     cos_phi .* sin_lambda,    -sin_phi .* sin_lambda,     cos_lambda;
    %     sin_phi,                   cos_phi,                   0
    % ];

    % This is the vectorized version of the transform matrix
    T = [
        transpose(cos_phi .* cos_lambda);
        transpose(cos_phi .* sin_lambda);
        transpose(sin_phi);
        transpose(-sin_phi .* cos_lambda);
        transpose(-sin_phi .* sin_lambda);
        transpose(cos_phi);
        transpose(-sin_lambda);
        transpose(cos_lambda);
        zeros(size(transpose(cos_lambda)));
    ];
    T = reshape(T, 3, 3, []);

    acceleration_spherical = transpose([a_r, a_phi, a_lambda]);
    acceleration_itrf = pagemtimes(T, permute(acceleration_spherical, [1, 3, 2]));
end

function [C_unnormalized, S_unnormalized] = UnnormalizeCoefficients(C, S)

    C_unnormalized = zeros(size(C));
    S_unnormalized = zeros(size(S));

    degree = size(S, 1) - 1;

    for ell = 0:1:degree
        % The coefficients are normalized, so the normalization factor must be applied to the legendre polynomials as
        % well.
        m_ = reshape(0:1:ell, 1, []);
        normalization_factor = sqrt(...
            (2 - (m_ == 0)) .* ...
            (2 .* ell + 1) .* ...
            factorial(ell - m_) ./ ...
            factorial(ell + m_));

        C_unnormalized(ell + 1, :) = C(ell + 1, :) .* [normalization_factor, zeros(1, degree + 1 - numel(normalization_factor))];
        S_unnormalized(ell + 1, :) = S(ell + 1, :) .* [normalization_factor, zeros(1, degree + 1 - numel(normalization_factor))];
    end
end