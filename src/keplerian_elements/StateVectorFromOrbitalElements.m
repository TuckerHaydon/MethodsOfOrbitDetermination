function [...
        position_eci, ...
        velocity_eci] = ...
    StateVectorFromOrbitalElements( ...
        semiparameter, ...
        eccentricity, ...
        inclination, ...
        longitude_of_ascending_node, ...
        argument_of_periapsis, ...
        true_anomaly, ...
        argument_of_latitude, ...
        true_longitude, ...
        true_longitude_of_periapsis, ...
        gravitational_parameter)
    % Converts the the orbital elements into ECI state vector position/velocity values.
    %
    % Requires:
    % - semiparameter:
    %   - (1, 1) double.
    %   - The orbit semiparameter
    %   - Unit meters.
    % - eccentricity:
    %   - (1, 1) double.
    %   - The orbit eccentricity.
    %   - Unitless.
    % - inclination:
    %   - (1, 1) double.
    %   - The orbit inclination.
    %   - Unit radians.
    % - longitude_of_ascending_node:
    %   - (1, 1) double.
    %   - The orbit longitude of the ascending node.
    %   - Unit radians.
    %   - This is set to NaN if it is undefined (equatorial orbit).
    % - argument_of_periapsis:
    %   - (1, 1) double.
    %   - The orbit argument of periapsis.
    %   - Unit radians.
    %   - This is set to NaN if it is undefined (circular orbit).
    % - true_anomaly:
    %   - (1, 1) double.
    %   - The orbit true anomaly.
    %   - Unit radians.
    %   - This is set to NaN if it is undefined (circular orbit).
    % - argument_of_latitude:
    %   - (1, 1) double.
    %   - The orbit argument of latitude.
    %   - Unit radians.
    %   - This is set to NaN if it is undefined (equatorial orbit).
    % - true_longitude:
    %   - (1, 1) double.
    %   - The orbit true longitude.
    %   - Unit radians.
    % - true_longitude_of_periapsis:
    %   - (1, 1) double.
    %   - This orbit true longitude of periapsis.
    %   - Unit radians.
    % - gravitational_parameter:
    %   - (1, 1) double.
    %   - The gravitational parameter of the primary body.
    %   - Units meters^3 / seconds^2.
    % 
    % Returns:
    % - position_eci:
    %   - (3, 1) double.
    %   - The position of the satellite in the ECI frame at epoch.
    %   - Unit meters.
    % - velocity_eci:
    %   - (3, 1) double.
    %   - The velocity of the satellite in the ECI frame at epoch.
    %   - Unit meters / second.
    %
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications
    arguments(Input)
        semiparameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
        eccentricity(1, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        inclination(1, 1) double {mustBeReal, mustBeFinite, mustBeInZeroToPi}
        longitude_of_ascending_node(1, 1) double {mustBeReal, mustBeNonnegativeOrNaN}
        argument_of_periapsis(1, 1) double {mustBeReal, mustBeNonnegativeOrNaN}
        true_anomaly(1, 1) double {mustBeReal, mustBeNonnegativeOrNaN}
        argument_of_latitude(1, 1) double {mustBeReal, mustBeNonnegativeOrNaN}
        true_longitude(1, 1) double {mustBeReal, mustBeFinite, mustBeInZeroTo2Pi}
        true_longitude_of_periapsis(1, 1) double {mustBeReal, mustBeFinite, mustBeInZeroTo2Pi}
        gravitational_parameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    arguments(Output)
        position_eci(3, 1) double
        velocity_eci(3, 1) double
    end

    % Rename
    p           = semiparameter;
    e_mag       = eccentricity;
    I           = inclination;
    Omega       = longitude_of_ascending_node;
    omega       = argument_of_periapsis;
    nu          = true_anomaly;
    u           = argument_of_latitude;
    lambda      = true_longitude;
    omega_tilde = true_longitude_of_periapsis;
    mu          = gravitational_parameter;

    % Quick analysis
    orbit_is_circular   = IsCircularOrbit(eccentricity);
    orbit_is_elliptical = IsEllipticalOrbit(eccentricity);
    orbit_is_equatorial = IsEquatorialOrbit(inclination);

    if orbit_is_circular
        % If the orbit is circular and equatorial, set the longitude of the ascending node and the argument of periapsis
        % to zero. Also set the true anomaly to the true longitude.
        if orbit_is_equatorial
            omega = 0;
            Omega = 0;
            nu    = lambda;
        % If the orbit is circular and inclined, set the argument of periapsis to 0. Also set the true anomaly to the
        % argument of latitude.
        else
            omega = 0;
            nu    = u;
        end
    elseif orbit_is_elliptical && orbit_is_equatorial
        % If the orbit is elliptical and equatorial, set the longitude of the ascending node to 0 and the argument of
        % periapsis to the true argument of periapsis.
        Omega = 0;
        omega = omega_tilde;
    end

    % Precompute quantities. Recommended by Vallado.
    cos_nu      = cos(nu);
    sin_nu      = sin(nu);
    mu_over_p   = mu / p;

    % Compute the position/velocity in the perifocal frame.
    % Vallado, Algorithm 10
    r_vec_pqw = [
        p .* cos_nu ./ (1 + e_mag .* cos_nu);
        p .* sin_nu ./ (1 + e_mag .* cos_nu);
        0
    ];

    v_vec_pqw = [
        -sqrt(mu_over_p) .* sin_nu;
        +sqrt(mu_over_p) .* (e_mag + cos_nu);
        0;
    ];

    C_pqw2ijk = ROT3(-Omega) * ROT1(-I) * ROT3(-omega);

    r_vec_ijk = C_pqw2ijk * r_vec_pqw;
    v_vec_ijk = C_pqw2ijk * v_vec_pqw;
    
    % Output
    position_eci = r_vec_ijk;
    velocity_eci = v_vec_ijk;
end

%% Helper functions
function [R3] = ROT3(alpha)
    % Defined by Vallado 3-15
    arguments(Input)
        alpha(1, 1) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        R3(3, 3) double {mustBeReal, mustBeFinite}
    end

    R3 = [
        +cos(alpha),   +sin(alpha),     0;
        -sin(alpha),   +cos(alpha),     0;
         0,             0,              1;
    ];
end

function [R1] = ROT1(alpha)
    % Defined by Vallado 3-15
    arguments(Input)
        alpha(1, 1) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        R1(3, 3) double {mustBeReal, mustBeFinite}
    end

    R1 = [
        1,   0,              0;
        0,  +cos(alpha),    +sin(alpha);
        0,  -sin(alpha),    +cos(alpha);
    ];
end