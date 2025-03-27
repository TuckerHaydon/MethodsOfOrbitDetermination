function [gravitational_acceleration_gcrf] = ...
    SimpleTwoBodyOffsetGravitationalAccelerationGCRF(...
        offset_gcrf, ...
        gravitational_parameter)
    % Computes the gravitational acceleration in the GCRF frame due to a non-Earth body (Sun or Moon).
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
    % - gravitational_acceleration_gcrf:
    %   - (3, :) double.
    %   - The gravitational acceleration vectors in the GCRF frame.
    %   - Unit m/s^2.
    arguments(Input)
        offset_gcrf(3, :) double {mustBeNumeric, mustBeReal}
        gravitational_parameter(1, 1) double {mustBeNumeric, mustBeReal}
    end

    arguments(Output)
        gravitational_acceleration_gcrf(3, :) double {mustBeNumeric, mustBeReal}
    end

    % Rename
    r_vec = offset_gcrf;
    mu    = gravitational_parameter;

    % Compute the two-body gravitation.
    r_mag = vecnorm(r_vec, 2, 1);

    gravitational_acceleration_gcrf = -mu .* r_vec ./ (r_mag.^3);

    % Output check
    assert(all(size(offset_gcrf) == size(gravitational_acceleration_gcrf)));
end