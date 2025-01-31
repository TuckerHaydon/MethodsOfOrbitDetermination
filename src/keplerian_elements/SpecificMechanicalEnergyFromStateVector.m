function [specific_mechanical_energy] = SpecificMechanicalEnergyFromStateVector(...
        position_eci, ...
        velocity_eci, ...
        gravitational_parameter)
    % Converts the state vector position/velocity values into the six orbital elements.
    %
    % Requires:
    % - position_eci:
    %   - (3, :) double.
    %   - The position of the satellite in the ECI frame.
    %   - Unit meters.
    % - velocity_eci:
    %   - (3, :) double.
    %   - The velocity of the satellite in the ECI frame.
    %   - Unit meters / second.
    % - gravitational_parameter:
    %   - (1, 1) double.
    %   - The gravitation parameter of the primary body.
    %   - Units m^3 / s^2.
    % 
    % Returns:
    % - specific_mechanical_energy:
    %   - (:, 1) double.
    %   - The orbit specific mechanical energy.
    %   - Unit m^2 / s^4.
    %   - Note, this can be negative.
    %
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications
    arguments(Input)
        position_eci(3, :) double {mustBeReal, mustBeFinite}
        velocity_eci(3, :) double {mustBeReal, mustBeFinite}
        gravitational_parameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    arguments(Output)
        specific_mechanical_energy(:, 1) double {mustBeReal, mustBeFinite}
    end

    % Rename
    r_vec = position_eci;
    v_vec = velocity_eci;
    mu = gravitational_parameter;

    % Vallado, 1-20
    xi = 0.5 .* reshape(dot(v_vec, v_vec, 1), [], 1) - mu ./ sqrt(reshape(dot(r_vec, r_vec, 1), [], 1));

    % Output
    specific_mechanical_energy = xi;
end