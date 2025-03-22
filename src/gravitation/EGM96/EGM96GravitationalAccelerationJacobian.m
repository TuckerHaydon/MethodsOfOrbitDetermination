function [gravitational_acceleration_jacobian_itrf] = EGM96GravitationalAccelerationJacobian(...
        position_itrf, ...
        earth_gravitational_parameter, ...
        earth_radius)
    % Computes the Jacobian of the ITRF J2 gravitational model
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
    % - gravitational_acceleration_jacobian_itrf:
    %   - (3, 3, :) double.
    %   - The Jacobian of the gravitational potential acceleration.
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
    end

    arguments(Output)
        gravitational_acceleration_jacobian_itrf(3, 3, :) double {mustBeReal, mustBeFinite}
    end

    % Load the table into memory
    persistent GM Re C S
    if isempty(GM) || isempty(Re) || isempty(C) || isempty(S)
        [GM, Re, C, S] = EGM96GravitationalModel20();
    end

    % Check that the provided coefficients are the same as the model coefficients.
    assert(abs(earth_gravitational_parameter - GM) < 1e-9, "Model mismatch!");
    assert(abs(earth_radius - Re) < 1e-9, "Model mismatch!");

    % The 'J2' constant is the C_bar(2, 0) value multiplied by -sqrt(5).
    J2 = -C(2 + 1, 0 + 1) * sqrt(5);

    % Rename
    mu = earth_gravitational_parameter;
    R0 = earth_radius;

    % From Haydon, Eqs 15 - 18
    u_hat_x = [1; 0; 0];
    u_hat_y = [0; 1; 0];
    u_hat_z = [0; 0; 1];

    % This should be a [3, N] array
    r     = vecnorm(position_itrf, 2, 1);
    r_hat = position_itrf ./ r;

    r_hat_x = r_hat(1, :);
    r_hat_y = r_hat(2, :);
    r_hat_z = r_hat(3, :);

    % Make vectors column vectors. These are [N, 1] vectors.
    r       = reshape(r, [], 1);
    r_hat_x = reshape(r_hat_x, [], 1);
    r_hat_y = reshape(r_hat_y, [], 1);
    r_hat_z = reshape(r_hat_z, [], 1);

    r_hat_z2 = r_hat_z.^2;

    % Precompute constants.
    % These are [N, 1] vectors
    mu_over_r3 = mu ./ r.^3;
    C1         = (3/2) .* J2 .* (R0 ./ r).^2;

    % This should be a [3, N] array
    dx_dr = transpose(mu_over_r3) .* (3 .* transpose(r_hat_x) .* r_hat - u_hat_x + transpose(C1) .* (...
        5 .* transpose(r_hat_x - 7 .* r_hat_z2 .* r_hat_x) .* r_hat + ...
        transpose(5 .* r_hat_z2 - 1) .* u_hat_x + ...
        10 .* transpose(r_hat_z .* r_hat_x) .* u_hat_z));

    dy_dr = transpose(mu_over_r3) .* (3 .* transpose(r_hat_y) .* r_hat - u_hat_y + transpose(C1) .* (...
        5 .* transpose(r_hat_y - 7 .* r_hat_z2 .* r_hat_y) .* r_hat + ...
        transpose(5 .* r_hat_z2 - 1) .* u_hat_y + ...
        10 .* transpose(r_hat_z .* r_hat_y) .* u_hat_z));

    dz_dr = transpose(mu_over_r3) .* (3 .* transpose(r_hat_z) .* r_hat - u_hat_z + transpose(C1) .* (...
        5 .* transpose(r_hat_z - 7 .* r_hat_z2 .* r_hat_z) .* r_hat + ...
        transpose(5 .* r_hat_z2 - 3) .* u_hat_z + ...
        10 .* transpose(r_hat_z .* r_hat_z) .* u_hat_z));

    G = [dx_dr; dy_dr; dz_dr];
    G = reshape(G, 3, 3, []);

    % Rename
    gravitational_acceleration_jacobian_itrf = G;
end