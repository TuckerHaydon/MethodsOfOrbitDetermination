function [gravitation_ecef] = ECEFJ2Gravitation( ...
        position_ecef, ...
        gravitational_parameter, ...
        dynamical_form_factor, ...
        semimajor_axis)
    % Computes the gravitational acceleration vector in the ECEF frame for satellite orbiting an ellipsoid.
    % 
    % Requires:
    % - position_ecef:
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
    % - gravitation_ecef:
    %   - (3, :) double.
    %   - The gravitational acceleration vector.
    % 
    % References:
    % - https://ieeexplore.ieee.org/document/7081494?arnumber=7081494
    arguments(Input)
        position_ecef(3, :) double {mustBeReal, mustBeFinite}
        gravitational_parameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
        dynamical_form_factor(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
        semimajor_axis(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    % Additional checks at end
    arguments(Output)
        gravitation_ecef(3, :) double {mustBeReal, mustBeFinite}
    end

    % Same function as ECI, just with vector changed.
    [gravitation_ecef] = ECIJ2Gravitation( ...
        position_ecef, ...
        gravitational_parameter, ...
        dynamical_form_factor, ...
        semimajor_axis);
end