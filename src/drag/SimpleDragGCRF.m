function [drag_gcrf] = SimpleDragGCRF(...
        position_gcrf, ...
        velocity_gcrf, ...
        coefficient_of_drag, ...
        effective_area, ...
        mass, ...
        reference_air_density, ...
        reference_radial_distance, ...
        decay_rate, ...
        earth_rotation_rate)
    % Produces a simple exponential drag model.
    %
    % The force of drag is modeled as:
    %   
    %   f = -(1/2) C_{D} (A_{eff} / m) rho_{A} |V_{A}| V_{A}
    %
    % where f is the force of drag in an inertial frame (the geocentric celestial reference frame; GCRF), C_{D} is the
    % coefficient of drag, A_{eff} is the effective area subject to the drag, m is the mass of the object experiencing
    % drag, rho_{A} is the atmospheric density, and V_{A} is the velocity vector of the object relative to the
    % atmosphere.
    %
    % The atmospheric density is a function of the radial distance of the satellite from the center of the Earth:
    % 
    %   rho_{A} = rho_{0} exp[-(r - r_{0}) / H]
    % 
    % Here, r_{0} is the reference altitude, rho_{0} is the atmospheric density at the reference altitude r_{0}, r is
    % the radial distance the object is from the center of the Earth, and H is a scale factor.
    %
    % Requires:
    % - position_gcrf:
    %   - (3, :) double.
    %   - The position of the satellite in the GCRF frame.
    %   - Unit m.
    % - velocity_gcrf:
    %   - (3, :) double.
    %   - The velocity of the satellite in the GCRF frame.
    %   - Unit m/s.
    % - coefficient_of_drag:
    %   - (1, 1) double.
    %   - The coefficient of drag.
    %   - Unitless.
    % - effective_area:
    %   - (1, 1) double.
    %   - The effective area of the satellite subject to drag.
    %   - Unit m^2.
    % - mass:
    %   - (1, 1) double.
    %   - The mass of the satellite.
    %   - Unit kg.
    % - reference_air_density:
    %   - (1, 1) double.
    %   - The reference air density.
    %   - Unit kg / m^3
    % - reference_radial_distance:
    %   - (1, 1) double.
    %   - The radial distance associated with the reference altitude.
    %   - Unit m.
    % - decay_rate:
    %   - (1, 1) double.
    %   - The atmospheric density decay rate.
    %   - Unit m.
    % - earth_rotation_rate:
    %   - (1, 1) double.
    %   - The instantaneous Earth rotation rate.
    %   - Unit rad / s.
    % 
    % Returns:
    % - drag_gcrf:
    %   - (3, :) double.
    %   - The drag expressed in the GCRF frame.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th Edition
    arguments(Input)
        position_gcrf(3, :) double {mustBeReal, mustBeFinite}
        velocity_gcrf(3, :) double {mustBeReal, mustBeFinite}
        coefficient_of_drag(:, 1) double {mustBeReal, mustBeFinite}
        effective_area(1, 1) double {mustBeReal, mustBeFinite}
        mass(1, 1) double {mustBeReal, mustBeFinite}
        reference_air_density(1, 1) double {mustBeReal, mustBeFinite}
        reference_radial_distance(1, 1) double {mustBeReal, mustBeFinite}
        decay_rate(1, 1) double {mustBeReal, mustBeFinite}
        earth_rotation_rate(1, 1) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        drag_gcrf(3, :) double {mustBeReal, mustBeFinite}
    end

    num_vectors = size(position_gcrf, 2);

    assert(all(size(velocity_gcrf) == size(position_gcrf)), "Vectors must be the same size!");
    assert(numel(coefficient_of_drag) == num_vectors);

    % Rename
    r_vec     = position_gcrf;
    v_vec     = velocity_gcrf;
    C_d       = coefficient_of_drag;
    A_eff     = effective_area;
    m         = mass;
    rho_0     = reference_air_density;
    r_0       = reference_radial_distance;
    H         = decay_rate;
    theta_dot = earth_rotation_rate;

    % Useful quantities
    % Size [N, 1]
    r_mag = vecnorm(r_vec, 2, 1);
    r_mag = reshape(r_mag, [], 1);

    % Compute the atmospheric density
    % Size [N, 1]
    rho_A = rho_0 .* exp(-(r_mag - r_0) ./ H);

    % Compute the velocity of the satellite relative to the atmopshere. Here, we're assuming that the atmosphere is
    % fixed to the surface of the Earth and that the GCRF and Earth-fixed frames are related via a simple Z-rotation.
    % Size [3, N]
    V_A_vec = v_vec + [
        theta_dot .* r_vec(2, :); 
       -theta_dot .* r_vec(1, :); 
        zeros(size(r_vec(1, :)));
    ];

    % Size [N, 1]
    V_A_mag = reshape(vecnorm(V_A_vec, 2, 1), [], 1);

    % Compute the drag
    % Size [3, N]
    drag_gcrf = -1/2 .* transpose(C_d) .* (A_eff ./ m) .* transpose(rho_A) .* transpose(V_A_mag) .* V_A_vec;
end