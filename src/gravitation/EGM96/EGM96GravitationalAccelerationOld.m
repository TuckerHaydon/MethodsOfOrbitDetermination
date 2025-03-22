function [acceleration_itrf] = EGM96GravitationalAccelerationOld(...
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

    mu_over_r  = mu ./ r;
    mu_over_r2 = mu ./ (r.^2);
    Re_over_r  = Re ./ r;

    % Geocentric longitude
    lambda = atan2(position_itrf(2, :), position_itrf(1, :));
    lambda = reshape(lambda, [], 1);

    % Geocentric latitude
    phi = asin(reshape(position_itrf(3, :), [], 1) ./ r);

    sin_phi = sin(phi);
    tan_phi = tan(phi);

    % Precompute the cos(m * lambda) and sin(m * lambda) terms.
    m = reshape(0:1:degree, 1, []);
    cos_m_lambda = cos(m .* lambda);
    sin_m_lambda = sin(m .* lambda);

    partial_R_partial_r      = zeros(num_vectors, 1);
    partial_R_partial_phi    = zeros(num_vectors, 1);
    partial_R_partial_lambda = zeros(num_vectors, 1);

    % Loops are Vallado 8-25
    for ell = 2:1:degree

        % Note that the matlab legendre function is off by a factor of ((-1).^m_) from its definition in gravitational
        % spherical harmonic literature.
        m_ = reshape(0:1:ell, 1, []);
        P = transpose(legendre(ell, sin_phi)) .* (-1).^m_;

        partial_R_partial_r_row      = zeros(num_vectors, 1);
        partial_R_partial_phi_row    = zeros(num_vectors, 1);
        partial_R_partial_lambda_row = zeros(num_vectors, 1);

        for m = 0:1:ell

            C_ell_m = C_unnormalized(ell + 1, m + 1);
            S_ell_m = S_unnormalized(ell + 1, m + 1);
            P_ell_m = P(m + 1);
            
            partial_R_partial_r_row = partial_R_partial_r_row + ...
                 P_ell_m .* (C_ell_m .* cos_m_lambda(:, m + 1) + S_ell_m .* sin_m_lambda(:, m + 1));

            if m + 2 <= ell
                P_ell_mp1 = P(m + 2);

                partial_R_partial_phi_row = partial_R_partial_phi_row + ...
                    (P_ell_mp1 - m .* tan_phi .* P_ell_m) .* ...
                    (C_ell_m .* cos_m_lambda(:, m + 1) + S_ell_m .* sin_m_lambda(:, m + 1));
            end

            partial_R_partial_lambda_row = partial_R_partial_lambda_row + ...
                m .* P_ell_m .* (S_ell_m .* cos_m_lambda(:, m + 1) - C_ell_m .* sin_m_lambda(:, m + 1));
        end

        Re_over_r_to_ell = power(Re_over_r, ell);

        partial_R_partial_r_row = Re_over_r_to_ell .* (ell + 1) .* partial_R_partial_r_row;
        partial_R_partial_r     = partial_R_partial_r            + partial_R_partial_r_row;

        partial_R_partial_phi_row = Re_over_r_to_ell     .* partial_R_partial_phi_row;
        partial_R_partial_phi     = partial_R_partial_phi + partial_R_partial_phi_row;

        partial_R_partial_lambda_row = Re_over_r_to_ell        .* partial_R_partial_lambda_row;
        partial_R_partial_lambda     = partial_R_partial_lambda + partial_R_partial_lambda_row;
    end

    partial_R_partial_r      = -mu_over_r2 .* partial_R_partial_r;
    partial_R_partial_phi    = mu_over_r   .* partial_R_partial_phi;
    partial_R_partial_lambda = mu_over_r   .* partial_R_partial_lambda;

    % Vallado, 8-27
    rI = position_itrf(1, :);
    rJ = position_itrf(2, :);
    rK = position_itrf(3, :);

    aI = ((1 ./ r) .* partial_R_partial_r - rK ./ (r.^2 .* sqrt(rI.^2 + rJ.^2)) .* partial_R_partial_phi) .* rI ...
        - (1 / (rI.^2 + rJ.^2) .* partial_R_partial_lambda) .* rJ ;

    aJ = ((1 ./ r) .* partial_R_partial_r - rK ./ (r.^2 .* sqrt(rI.^2 + rJ.^2)) .* partial_R_partial_phi) .* rJ ...
        + (1 ./ (rI.^2 + rJ.^2) .* partial_R_partial_lambda) .* rI;

    aK = (1 ./ r) .* partial_R_partial_r .* rK + (sqrt(rI.^2 + rJ.^2) ./ r.^2) .* partial_R_partial_phi;

    acceleration_itrf = [aI; aJ; aK];

    % Add two-body
    acceleration_itrf = acceleration_itrf - mu .* position_itrf ./ r.^3;
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