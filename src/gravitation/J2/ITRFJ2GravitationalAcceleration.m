function [gravitational_acceleration_itrf] = ITRFJ2GravitationalAcceleration( ...
        position_itrf, ...
        earth_gravitational_parameter, ...
        earth_radius, ...
        J2)
    % Computes the Earth gravitational acceleration vector using a J2 model.
    % 
    % Requires:
    % - position_itrf:
    %   - (3, :) double.
    %   - The position at which to calculate the gravitational acceleration vector in the ITRF frame.
    %   - Unit m.
    % - earth_gravitational_parameter:
    %   - (1, 1) double.
    %   - The Earth gravitational parameter.
    %   - Units meters^3 / seconds^2.
    %   - Unitless.
    % - earth_radius:
    %   - (1, 1) double.
    %   - The Earth radius.
    %   - Unit meters.
    % - J2:
    %   - (1, 1) double.
    %   - The J2 dynamical form factor.
    % 
    % Returns:
    % - gravitational_acceleration_itrf:
    %   - (3, :) double.
    %   - The gravitational acceleration vector in the ITRF frame.
    % 
    % References:
    % - https://ieeexplore.ieee.org/document/7081494?arnumber=7081494
    arguments(Input)
        position_itrf(3, :) double {mustBeReal, mustBeFinite}
        earth_gravitational_parameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
        earth_radius(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
        J2(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    % Additional checks at end
    arguments(Output)
        gravitational_acceleration_itrf(3, :) double {mustBeReal, mustBeFinite}
    end

    % Rename
    r_vec = position_itrf;
    mu    = earth_gravitational_parameter;
    Re    = earth_radius;

    % Compute the two-body gravitation.
    r_mag = vecnorm(r_vec, 2, 1);

    gravitational_acceleration_itrf = ...
        -(mu  ./ (r_mag.^3)) .* (r_vec + (3/2) .* J2 .* (Re ./ r_mag).^2 .* [
            (1 - 5 .* (r_vec(3,:) ./ r_mag).^2) .* r_vec(1, :);
            (1 - 5 .* (r_vec(3,:) ./ r_mag).^2) .* r_vec(2, :);
            (3 - 5 .* (r_vec(3,:) ./ r_mag).^2) .* r_vec(3, :);
    ]);

    % Output checks 
    assert(size(gravitational_acceleration_itrf, 2) == size(position_itrf, 2));
end