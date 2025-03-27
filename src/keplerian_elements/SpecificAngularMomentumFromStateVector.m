function [specific_angular_momentum] = SpecificAngularMomentumFromStateVector( ...
        position_eci, ...
        velocity_eci)
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
    % 
    % Returns:
    % - specific_angular_momentum:
    %   - (3, :) double.
    %   - The orbit specific angular momentum.
    %   - Unit m^2 / s.
    %
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications
    arguments(Input)
        position_eci(3, :) double {mustBeReal, mustBeFinite}
        velocity_eci(3, :) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        specific_angular_momentum(3, :) double {mustBeReal, mustBeFinite}
    end

    % Rename
    r_vec = position_eci;
    v_vec = velocity_eci;

    % Vallado, 1-15
    h_vec = cross(r_vec, v_vec, 1);

    % Output
    specific_angular_momentum = h_vec;

    % Output check
    assert(size(specific_angular_momentum, 2) == size(position_eci, 2));
end