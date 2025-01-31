function [angular_velocity_eci] = AngularVelocityFromStateVector( ...
        position_eci, ...
        velocity_eci)
    % Computes the instantaneous angular velocity vector from the ECI state vector.
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
    % - angular_velocity_eci:
    %   - (3, :) double.
    %   - The instantaneous angular velocity in rad/s.
    % 
    % References:
    % - https://phys.libretexts.org/Bookshelves/University_Physics/Book%3A_Introductory_Physics_-_Building_Models_to_Describe_Our_World_(Martin_Neary_Rinaldo_and_Woodman)/11%3A_Rotational_dynamics/11.01%3A_Rotational_kinematic_vectors
    arguments(Input)
        position_eci(3, :) double {mustBeReal, mustBeFinite}
        velocity_eci(3, :) double {mustBeReal, mustBeFinite}
    end

    % Additional checks at end
    arguments(Output)
        angular_velocity_eci(3, :) double {mustBeReal, mustBeFinite}
    end

    % Rename
    r = position_eci;
    v = velocity_eci;

    omega = cross(r, v, 1) ./ dot(r, r, 1);

    % Output
    angular_velocity_eci = omega;

    % Checks
    assert(size(angular_velocity_eci, 2) == size(position_eci, 2));
end