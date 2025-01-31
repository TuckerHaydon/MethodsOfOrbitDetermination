function [gravitation_eci] = ECISimpleTwoBodyGravitation( ...
        position_eci, ...
        gravitational_parameter)
    % - gravitational_parameter:
    %   - (1, 1) double.
    %   - The gravitational parameter of the primary body.
    %   - Units meters^3 / seconds^2.
    arguments(Input)
        position_eci(3, :) double {mustBeReal, mustBeFinite}
        gravitational_parameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    % Additional checks at end
    arguments(Output)
        gravitation_eci(3, :) double {mustBeReal, mustBeFinite}
    end

    % Rename
    r_vec = position_eci;
    mu    = gravitational_parameter;

    % Compute the two-body gravitation.
    r_mag = vecnorm(r_vec, 2, 1);

    gravitation_eci = -mu .* r_vec ./ (r_mag.^3);

    % Output checks 
    assert(size(gravitation_eci, 2) == size(position_eci, 2));
end