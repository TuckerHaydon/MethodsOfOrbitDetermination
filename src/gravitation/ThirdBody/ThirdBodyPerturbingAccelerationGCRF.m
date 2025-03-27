function [gravitational_acceleration_gcrf] = ...
    ThirdBodyPerturbingAccelerationGCRF(...
        satellite_position_gcrf, ...
        perturbing_body_position_gcrf, ...
        gravitational_parameter)
    % Computes the gravitational acceleration in the GCRF frame due to a non-Earth body (Sun or Moon).
    %
    % Requires:
    % - satellite_position_gcrf:
    %   - (3, :) double.
    %   - The position of the satellite in the GCRF frame.
    %   - Unit meters
    % - perturbing_body_position_gcrf:
    %   - (3, :) double.
    %   - The position of the perturbing body in the GCRF frame.
    %   - Unit meters
    % - gravitational_parameter:
    %   - (1, 1) double.
    %   - The gravitational parameter of the non-Earth body.
    %   - Unit m^3 / s^2.
    %
    % Returns:
    % - gravitational_acceleration_gcrf:
    %   - (3, :) double.
    %   - The perturbing gravitational acceleration vectors in the GCRF frame.
    %   - Unit m/s^2.
    %
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th edition
    arguments(Input)
        satellite_position_gcrf(3, :) double {mustBeNumeric, mustBeReal}
        perturbing_body_position_gcrf(3, 1) double {mustBeNumeric, mustBeReal}
        gravitational_parameter(1, 1) double {mustBeNumeric, mustBeReal}
    end

    arguments(Output)
        gravitational_acceleration_gcrf(3, :) double {mustBeNumeric, mustBeReal}
    end

    % Rename
    r_vec_sat  = satellite_position_gcrf;
    r_vec_body = perturbing_body_position_gcrf;
    mu         = gravitational_parameter;

    % Need the vector from the body to the satellite/Earth
    r_vec_body2sat = r_vec_sat - r_vec_body;
    r_vec_body2earth = -r_vec_body;

    % Compute the two-body gravitation.
    r_mag_body2sat  = vecnorm(r_vec_body2sat, 2, 1);
    r_mag_body2earth = vecnorm(r_vec_body2earth, 2, 1);

    % Vallado, 8-33
    gravitational_acceleration_gcrf = -mu .* (r_vec_body2sat ./ r_mag_body2sat.^3 - r_vec_body2earth ./ r_mag_body2earth.^3);

    % Output check
    assert(all(size(satellite_position_gcrf) == size(gravitational_acceleration_gcrf)));
end