function [H] = ApparentRangeAndRangeRateJacobian( ...
        satellite_position_at_rx_gcrf, ...
        satellite_velocity_at_rx_gcrf, ...
        station_position_at_rx_gcrf, ...
        station_velocity_at_rx_gcrf, ...
        light_time)
    % Computes the Jacobian of apparent range and range-rate measurements with respect to the satellite's GCRF position
    % at the time of measurement observation.
    %
    % Requires:
    % - satellite_position_at_rx_gcrf:
    %   - (3, 1) double.
    %   - The satellite position in the GCRF frame at the time of observation.
    %   - Unit meters.
    % - satellite_velocity_at_rx_gcrf:
    %   - (3, 1) double.
    %   - The satellite velocity in the GCRF frame at the time of observation.
    %   - Unit meters / second.
    % - station_position_at_rx_gcrf:
    %   - (3, 1) double.
    %   - The station position in the GCRF frame at the time of observation.
    %   - Unit meters.
    % - station_velocity_at_rx_gcrf: 
    %   - (3, 1) double.
    %   - The station velocity in the GCRF frame at the time of observation.
    %   - Unit meters / second.
    % - light_time:
    %   - (1, 1) double.
    %   - The computed light time associated with the measurement.
    %   - Unit seconds.
    %
    % Returns:
    % - H:
    %   - (2, 6) double.
    %   - The apparent range and range-rate Jacobian w.r.t. satellite position and velocity.
    %
    % References:
    % - 
    arguments(Input)
        satellite_position_at_rx_gcrf(3, 1) double {mustBeReal, mustBeFinite}
        satellite_velocity_at_rx_gcrf(3, 1) double {mustBeReal, mustBeFinite}
        station_position_at_rx_gcrf(3, 1) double {mustBeReal, mustBeFinite}
        station_velocity_at_rx_gcrf(3, 1) double {mustBeReal, mustBeFinite}
        light_time(1, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative} 
    end

    arguments(Output)
        H(2, 6) double {mustBeReal, mustBeFinite}
    end

    satellite_position_at_tx_gcrf = satellite_position_at_rx_gcrf - light_time * satellite_velocity_at_rx_gcrf;
    satellite_velocity_at_tx_gcrf = satellite_velocity_at_rx_gcrf; % Constant velocity
    
    delta_position_vector = station_position_at_rx_gcrf - satellite_position_at_tx_gcrf;
    delta_velocity_vector = station_velocity_at_rx_gcrf - satellite_velocity_at_tx_gcrf;

    r_mag = norm(delta_position_vector);
    r_hat = delta_position_vector ./ r_mag;
    
    % Jacobian
    H(1, 1:3) = -r_hat;
    H(1, 4:6) = zeros(1, 3);

    H(2, 1:3) = -(eye(3) - r_hat * transpose(r_hat)) * delta_velocity_vector / r_mag;
    H(2, 4:6) = -r_hat;
end