function [gravitation_eci] = ECIJ2GravitationalAcceleration( ...
        position_eci, ...
        gravitational_parameter, ...
        dynamical_form_factor, ...
        semimajor_axis)
    % Computes the gravitational acceleration vector in the ECI frame for satellite orbiting an ellipsoid.
    % 
    % Requires:
    % - position_eci:
    %   - (3, :) double.
    %   - The position at which to calculate the gravitational acceleration vector.
    %   - Unit m.
    % - gravitational_parameter:
    %   - (1, 1) double.
    %   - The gravitational parameter of the primary body.
    %   - Units meters^3 / seconds^2.
    % - dynamical_form_factor:
    %   - (1, 1) double.
    %   - The J2 dynamical form factor.
    %   - Unitless.
    % - semimajor_axis:
    %   - (1, 1) double.
    %   - The ellipsoid semimajor axis.
    %   - Unit meters.
    % 
    % Returns:
    % - gravitation_eci:
    %   - (3, :) double.
    %   - The gravitational acceleration vector.
    % 
    % References:
    % - https://ieeexplore.ieee.org/document/7081494?arnumber=7081494
    arguments(Input)
        position_eci(3, :) double {mustBeReal, mustBeFinite}
        gravitational_parameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
        dynamical_form_factor(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
        semimajor_axis(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    % Additional checks at end
    arguments(Output)
        gravitation_eci(3, :) double {mustBeReal, mustBeFinite}
    end

    % Rename
    r_vec = position_eci;
    mu    = gravitational_parameter;
    J2    = dynamical_form_factor;
    Re    = semimajor_axis;

    % Compute the two-body gravitation.
    r_mag = vecnorm(r_vec, 2, 1);

    gravitation_eci = -(mu  ./ (r_mag.^3)) .* (r_vec + (3/2) .* J2 .* (Re ./ r_mag).^2 .* [
        (1 - 5 .* (r_vec(3,:) ./ r_mag).^2) .* r_vec(1, :);
        (1 - 5 .* (r_vec(3,:) ./ r_mag).^2) .* r_vec(2, :);
        (3 - 5 .* (r_vec(3,:) ./ r_mag).^2) .* r_vec(3, :);
    ]);

    % Output checks 
    assert(size(gravitation_eci, 2) == size(position_eci, 2));
end