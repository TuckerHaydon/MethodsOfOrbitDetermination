function [position_jacobian_gcrf] = ...
    SimpleTwoBodyOffsetGravitationalAccelerationJacobianGCRF(...
        offset_gcrf, ...
        gravitational_parameter)
    % Computes the Jacobian of the gravitational acceleration vectors due to offset bodies (Sun or Moon) in the GCRF frame.
    %
    % Requires:
    % - offset_gcrf:
    %   - (3, :) double.
    %   - The position vector from the object to the satellite in the GCRF frame.
    %   - Unit meters
    % - gravitational_parameter:
    %   - (1, 1) double.
    %   - The gravitational parameter of the non-Earth body.
    %   - Unit m^3 / s^2.
    %
    % Returns:
    % - position_jacobian_gcrf:
    %   - (3, 3, :) double.
    %   - The position Jacobians of the gravitational acceleration vectors in the GCRF frame.
    %   - Unit (m/s^2) / m.
    arguments(Input)
        offset_gcrf(3, :) double {mustBeNumeric, mustBeReal}
        gravitational_parameter(1, 1) double {mustBeNumeric, mustBeReal}
    end

    arguments(Output)
        position_jacobian_gcrf(3, 3, :) double {mustBeNumeric, mustBeReal}
    end

    % Rename
    mu          = gravitational_parameter;
    delta_r_vec = offset_gcrf;

    % Compute the two-body gravitation.
    delta_r_mag = vecnorm(delta_r_vec, 2, 1);
    delta_r_hat = delta_r_vec ./ delta_r_mag;

    % Sizes:
    % delta_r_vec [3, N]
    % delta_r_mag [1, N]
    % delta_r_hat [3, N]
    position_jacobian_gcrf = -mu ./ permute(delta_r_mag, [1, 3, 2]).^3 .* (eye(3) - 3 .* permute(delta_r_hat, [1, 3, 2]) .* permute(delta_r_hat, [3, 1, 2])); 

    % Output check
    num_vectors = size(offset_gcrf, 2);
    assert(all(size(position_jacobian_gcrf, 3) == num_vectors));
end