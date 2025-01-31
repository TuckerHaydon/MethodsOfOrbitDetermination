function [...
        semiparameter, ...
        semimajor_axis, ...
        eccentricity, ...
        inclination, ...
        longitude_of_ascending_node, ...
        argument_of_periapsis, ...
        true_anomaly, ...
        argument_of_latitude, ...
        true_longitude, ...
        true_longitude_of_periapsis] = ...
    OrbitalElementsFromStateVector( ...
        position_eci, ...
        velocity_eci, ...
        gravitational_parameter)
    % Converts the state vector position/velocity values into the Keplerian orbital elements.
    %
    % Note, this function is not vectorized.
    %
    % Requires:
    % - position_eci:
    %   - (3, 1) double.
    %   - The position of the satellite in the ECI frame at epoch.
    %   - Unit meters.
    % - velocity_eci:
    %   - (3, 1) double.
    %   - The velocity of the satellite in the ECI frame at epoch.
    %   - Unit meters / second.
    % - gravitational_parameter:
    %   - (1, 1) double.
    %   - The gravitational parameter of the primary body.
    %   - Units meters^3 / seconds^2.
    % 
    % Returns:
    % - semiparameter:
    %   - (1, 1) double.
    %   - The orbit semiparameter
    %   - Unit meters.
    % - semimajor_axis:
    %   - (1, 1) double.
    %   - The orbit semimajor axis.
    %   - Unit meters.
    %   - This is set to NaN if it is undefined (parabolic orbit).
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
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications
    arguments(Input)
        position_eci(3, 1) double {mustBeReal, mustBeFinite}
        velocity_eci(3, 1) double {mustBeReal, mustBeFinite}
        gravitational_parameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
    end

    arguments(Output)
        semiparameter(1, 1) double {mustBeReal, mustBeFinite, mustBePositive}
        semimajor_axis(1, 1) double {mustBeReal, mustBePositiveOrNaN}
        eccentricity(1, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        inclination(1, 1) double {mustBeReal, mustBeFinite, mustBeInZeroToPi}
        longitude_of_ascending_node(1, 1) double {mustBeReal, mustBeNonnegativeOrNaN}
        argument_of_periapsis(1, 1) double {mustBeReal, mustBeNonnegativeOrNaN}
        true_anomaly(1, 1) double {mustBeReal, mustBeNonnegativeOrNaN}
        argument_of_latitude(1, 1) double {mustBeReal, mustBeNonnegativeOrNaN}
        true_longitude(1, 1) double {mustBeReal, mustBeFinite, mustBeInZeroTo2Pi}
        true_longitude_of_periapsis(1, 1) double {mustBeReal, mustBeFinite, mustBeInZeroTo2Pi}
    end

    % Rename
    r_vec = position_eci;
    v_vec = velocity_eci;
    mu    = gravitational_parameter;

    r_mag = norm(r_vec);
    r_hat = r_vec ./ r_mag;

    v_mag = norm(v_vec);

    % First find specific angular momentum. It is defined for all orbits.
    % https://en.wikipedia.org/wiki/Specific_angular_momentum
    h_vec = SpecificAngularMomentumFromStateVector(r_vec, v_vec);
    h_mag = norm(h_vec);

    % Next, find the node vector. It is defined for all orbits.
    k_hat = [0; 0; 1];
    n_vec = cross(k_hat, h_vec);
    n_mag = norm(n_vec);
    n_hat = n_vec ./ n_mag;

    % Next, find the eccentricity/Laplace vector. It is defined for all orbits.
    % Vallado, 2-80
    e_vec = (1.0 / mu) * ((v_mag.^2 - mu ./ r_mag) .* r_vec - dot(r_vec, v_vec, 1) .* v_vec);
    e_mag = norm(e_vec);
    e_hat = e_vec ./ e_mag;

    % Next, find the specific mechanical energy. It is defined for all orbits.
    xi = SpecificMechanicalEnergyFromStateVector(r_vec, v_vec, mu);

    % Inclination
    % Vallado, 2-84
    I = acos(h_vec(3) ./ h_mag);

    % Determine the shape of the orbit.
    orbit_is_equatorial = IsEquatorialOrbit(I);
    orbit_is_circular   = IsCircularOrbit(e_mag);
    orbit_is_parabolic  = IsParabolicOrbit(e_mag);

    % Determine the semiparameter
    % Vallado, 1-19
    p = h_mag.^2 ./ mu;

    % Determine the semimajor axis
    % Vallado, 1-21
    if ~orbit_is_parabolic
        a = -mu ./ (2.0 .* xi);
    else
        a = NaN;
    end

    % If orbit is equatorial, the longitude of the ascending node is undefined.
    % Vallado, 2-86
    if orbit_is_equatorial
        Omega = NaN;
    else
        Omega = acos(n_vec(1) / n_mag);

        if n_vec(2) < 0
            Omega = 2.0 .* pi - Omega;
        end
    end

    % If orbit is circular, the argument of periapsis is undefined.
    % Vallado, 2-87
    if orbit_is_circular
        omega = NaN;
    else
        omega = acos(dot(n_hat, e_hat));

        if e_vec(3) < 0
            omega = 2.0 .* pi - omega;
        end
    end

    % If orbit is circular, the true anomaly is undefined.
    % Vallado, 2-88
    if orbit_is_circular
        nu = NaN;
    else
        tmp = dot(e_hat, r_hat);
        if abs(tmp) > 1
            warning("OrbitalElementsFromStateVector: abs(dot(e_hat, r_hat)) > 1 by %0.3e. Clamping to unity. " + ...
                "This may indicate numerical instabilities, or might simply be a floating point issue.", abs(tmp) - 1);
            tmp = tmp ./ abs(tmp);
        end

        nu = acos(tmp);

        if dot(r_vec, v_vec) < 0
            nu = 2.0 .* pi - nu;
        end
    end

    % If an orbit is equatorial, the argument of latitude is undefined.
    % Vallado, 2-91
    if orbit_is_equatorial
        u = NaN;
    else
        u = acos(dot(n_hat, r_hat));

        if r_vec(3) < 0
            u = 2.0 .* pi - u;
        end
    end

    % The true longitude is defined for all orbits.
    % Vallado, 2-94
    lambda = acos(r_hat(1));

    if r_vec(2) < 0
        lambda = 2.0 .* pi - lambda;
    end

    % The true longitude of periapsis is defined for all orbits.
    % Vallado, 2-89
    omega_tilde = acos(e_hat(1));

    if e_vec(2) < 0
        omega_tilde = 2.0 .* pi - omega_tilde;
    end

    % Output
    semiparameter               = p;
    semimajor_axis              = a;
    eccentricity                = e_mag;
    inclination                 = I;
    longitude_of_ascending_node = Omega;
    argument_of_periapsis       = omega;
    true_anomaly                = nu;
    argument_of_latitude        = u;
    true_longitude              = lambda;
    true_longitude_of_periapsis = omega_tilde;
end