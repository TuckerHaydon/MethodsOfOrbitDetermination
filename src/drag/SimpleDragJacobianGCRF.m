function [position_gcrf_jacobian, velocity_gcrf_jacobian, coefficient_of_drag_jacobian] = ...
    SimpleDragJacobianGCRF(...
        position_gcrf, ...
        velocity_gcrf, ...
        coefficient_of_drag, ...
        effective_area, ...
        mass, ...
        reference_air_density, ...
        reference_radial_distance, ...
        decay_rate, ...
        earth_rotation_rate)
    arguments(Input)
        position_gcrf(3, :) double {mustBeReal, mustBeFinite}
        velocity_gcrf(3, :) double {mustBeReal, mustBeFinite}
        coefficient_of_drag(:, 1) double {mustBeReal, mustBeFinite}
        effective_area(1, 1) double {mustBeReal, mustBeFinite}
        mass(1, 1) double {mustBeReal, mustBeFinite}
        reference_air_density(1, 1) double {mustBeReal, mustBeFinite}
        reference_radial_distance(1, 1) double {mustBeReal, mustBeFinite}
        decay_rate(1, 1) double {mustBeReal, mustBeFinite}
        earth_rotation_rate(1, 1) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        position_gcrf_jacobian(3, 3, :) double {mustBeReal, mustBeFinite}
        velocity_gcrf_jacobian(3, 3, :) double {mustBeReal, mustBeFinite}
        coefficient_of_drag_jacobian(3, 1, :) double {mustBeReal, mustBeFinite}
    end

    num_vectors = size(position_gcrf, 2);

    assert(all(size(velocity_gcrf) == size(position_gcrf)), "Vectors must be the same size!");
    assert(numel(coefficient_of_drag) == num_vectors);


    % Rename
    r_vec     = position_gcrf;
    v_vec     = velocity_gcrf;
    C_d       = coefficient_of_drag;
    A_eff     = effective_area;
    m         = mass;
    rho_0     = reference_air_density;
    r_0       = reference_radial_distance;
    H         = decay_rate;
    theta_dot = earth_rotation_rate;

    num_vectors = size(position_gcrf, 2);

    % Size [N, 1]
    r_mag = vecnorm(r_vec, 2, 1);
    r_mag = reshape(r_mag, [], 1);

    % Size [3, N]
    r_hat = r_vec ./ transpose(r_mag);

    % Compute the atmospheric density.
    % Size [N, 1]
    rho_A = rho_0 .* exp(-(r_mag - r_0) ./ H);

    % Compute the velocity of the satellite relative to the atmopshere. Here, we're assuming that the atmosphere is
    % fixed to the surface of the Earth and that the GCRF and Earth-fixed frames are related via a simple Z-rotation.
    % Size [3, N]
    V_A_vec = v_vec + [
        theta_dot .* r_vec(2, :); 
       -theta_dot .* r_vec(1, :); 
        zeros(size(r_vec(1, :)));
    ];

    % Size [N, 1]
    V_A_mag = reshape(vecnorm(V_A_vec, 2, 1), [], 1);

    % Size [3, N]
    V_A_hat = V_A_vec ./ transpose(V_A_mag);

    % Partial derivative of rho_A w.r.t. position
    % Size [N, 3]
    drhoA_dr = -rho_0 ./ H .* exp(-(r_mag - r_0) ./ H) .* transpose(r_hat);

    % Partial derivative of airspeed vector w.r.t. position
    % Size [3, 3, N]
    dVAvec_dr = [
        zeros(1, num_vectors);
        -theta_dot .* ones(1, num_vectors);
        zeros(1, num_vectors);
        +theta_dot .* ones(1, num_vectors);
        zeros(1, num_vectors);
        zeros(1, num_vectors);
        zeros(1, num_vectors);
        zeros(1, num_vectors);
        zeros(1, num_vectors);
    ];
    dVAvec_dr = reshape(dVAvec_dr, 3, 3, []);

    % Size [N, 3]
    dVAmag_dr = permute(pagemtimes(permute(V_A_hat, [3, 1, 2]), dVAvec_dr), [3, 2, 1]);

    % Size [3, 3, N]
    dVAvec_dv = repmat(eye(3, 3), 1, 1, num_vectors);

    % Size [N, 3]
    dVAmag_dv = transpose(V_A_hat);

    % Now compute jacobians
    % The Jacobian w.r.t. the coefficient of drag
    coefficient_of_drag_jacobian = (-1/2 .* (A_eff ./ m)) .* permute(rho_A, [3, 2, 1]) .* permute(V_A_mag, [3, 2, 1]) .* permute(V_A_vec, [1, 3, 2]);

    % This works for a single vector.
    % position_gcrf_jacobian = ...
    %     -1/2 .* C_d .* (A_eff ./ m) .* drhoA_dr .* V_A_vec   .* V_A_mag ...
    %     -1/2 .* C_d .* (A_eff ./ m) .* rho_A    .* dVAvec_dr .* V_A_mag ...
    %     -1/2 .* C_d .* (A_eff ./ m) .* rho_A    .* V_A_vec   .* dVAmag_dr;

    % Sizes to reference for permute:
    % - C_d       [N, 1]
    % - drhoA_dr  [N, 3]
    % - V_A_vec   [3, N]
    % - V_A_mag   [N, 1]
    % - rho_A     [N, 1]
    % - dVAvec_dr [3, 3, N]
    % - dVAmag_dr [N, 3]
    % - dVAvec_dv [3, 3, N]
    % - dVAmag_dv [N, 3]
    position_gcrf_jacobian = ...
        -1/2 .* permute(C_d, [3, 2, 1]) .* (A_eff ./ m) .* permute(drhoA_dr, [3, 2, 1]) .* permute(V_A_vec, [1, 3, 2]) .* permute(V_A_mag, [3, 2, 1]) ...
        -1/2 .* permute(C_d, [3, 2, 1]) .* (A_eff ./ m) .* permute(rho_A, [3, 2, 1]) .* dVAvec_dr .* permute(V_A_mag, [3, 2, 1]) ...
        -1/2 .* permute(C_d, [3, 2, 1]) .* (A_eff ./ m) .* permute(rho_A, [3, 2, 1]) .* permute(V_A_vec, [1, 3, 2]) .* permute(dVAmag_dr, [3, 2, 1]);


    velocity_gcrf_jacobian = ...
        -1/2 .* permute(C_d, [3, 2, 1]) .* (A_eff ./ m) .* permute(rho_A, [3, 2, 1]) .* dVAvec_dv .* permute(V_A_mag, [3, 2, 1]) ...
        -1/2 .* permute(C_d, [3, 2, 1]) .* (A_eff ./ m) .* permute(rho_A, [3, 2, 1]) .* permute(V_A_vec, [1, 3, 2]) .* permute(dVAmag_dv, [3, 2, 1]);
end