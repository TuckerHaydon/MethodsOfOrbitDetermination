function [gravitational_potential] = ECIJ2GravitationalPotential( ...
        position_eci, ...
        gravitational_parameter, ...
        dynamical_form_factor, ...
        semimajor_axis)
    % Computes the gravitational potential from an ECI position vector.
    % 
    % Requires:
    % - position_eci:
    %   - (3, :) double.
    %   - The position at which to calculate the gravitational potential.
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
    % - gravitational_potential:
    %   - (:, 1) double.
    %   - The gravitational potential.
    arguments(Input)
        position_eci(3, :) double {mustBeReal, mustBeFinite}
        gravitational_parameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
        dynamical_form_factor(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
        semimajor_axis(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    % Additional checks at end
    arguments(Output)
        gravitational_potential(:, 1) double {mustBeReal, mustBeFinite}
    end

    % Rename
    r_vec = position_eci;
    mu    = gravitational_parameter;
    J2    = dynamical_form_factor;
    Re    = semimajor_axis;

    % Compute the two-body gravitation.
    r_mag = vecnorm(r_vec, 2, 1);

    sin_phi = r_vec(3, :) ./ r_mag;
    gravitational_potential = mu ./ r_mag .* (1 - J2 .* (Re ./ r_mag).^2 .* (3/2 .* sin_phi.^2 - 1/2));
    gravitational_potential = transpose(gravitational_potential);

    % Output checks 
    assert(size(gravitational_potential, 1) == size(position_eci, 2));
end