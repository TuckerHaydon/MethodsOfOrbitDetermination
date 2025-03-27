function [apparent_range, apparent_range_rate, light_time] = ApparentRangeAndRangeRate( ...
        satellite_position_at_rx_gcrf, ...
        satellite_velocity_at_rx_gcrf, ...
        station_position_at_rx_gcrf, ...
        station_velocity_at_rx_gcrf, ...
        options) 
    % Computes the apparent one-way range and range rate of a satellite as viewed by a ground-fixed station.
    %
    % Let t be the time a station receives a radio signal from a satellite, p_{station}^{GCRF}[t] be the position of an
    % observation station at time t, and p_{satellite}^{GCRF}[t] be the position of a satellite at time t. The apparent
    % one-way range of the satellite is
    %
    %   r = | p_{station}^{GCRF}[t] - p_{satellite}^{GCRF}[t - dt] |
    % 
    % where dt is the light time (the time it takes a radio signal to travel from the satellite to the monitoring
    % station). Here, p_{satellite}^{GCRF}[t - dt] is the position of the satellite when it sent the radio
    % communication.
    %
    % The prior position of the satellite p_{satellite}^{GCRF}[t - dt] could be computed by backwards integrating the
    % dynamics of the satellite. This would be computationally expensive. Alternatively, if the light time dt is small,
    % a first-order constant-velocity approximation could be invoked. In such a case, the apparent range reduces to
    %
    %   r = | p_{station}^{GCRF}[t] - (p_{satellite}^{GCRF}[t] - v_{satellite}^{GCRF}[t] * dt) |
    %
    % where v_{satellite}^{GCRF}[t] is the assumed-constant velocity of the satellite at time [t].
    % 
    % It is this first-order constant-velocity approximation that is used. 
    % 
    % The apparent range-rate is the time derivative of this apparent range expression:
    % 
    %   d/dt r = dot[\hat{r}[t, dt], (v_{station}^{GCRF}[t] - v_{satellite}^{GCRF}[t])]
    %
    % Here, \hat{r} [t, dt] is the unit vector of 
    %   (p_{station}^{GCRF} - (p_{satellite}^{GCRF}[t] - v_{satellite}^{GCRF}[t] * dt))
    %
    % Requires:
    % - satellite_position_at_rx_gcrf:
    %   - (3, 1) double.
    %   - The satellite position at the time of observation (at time t) in the GCRF frame.
    %   - Unit meters.
    % - satellite_velocity_at_rx_gcrf:
    %   - (3, 1) double.
    %   - The satellite velocity at the time of observation (at time t) in the GCRF frame.
    %   - Unit meters / second.
    % - station_position_at_rx_gcrf:
    %   - (3, 1) double.
    %   - The station position at the time of observation (at time t) in the GCRF frame.
    %   - Unit meters.
    % - station_velocity_at_rx_gcrf:
    %   - (3, 1) double.
    %   - The station velocity at the time of observation (at time t) in the GCRF frame.
    %   - Unit meters / second.
    % - options (optional struct):
    %   - threshold:
    %     - (1, 1) double.
    %     - Iteration threshold of the light second computation.
    % %   - Unit seconds.
    %     - Defaults to 1/10 of a meter in light time.
    % 
    % Returns:
    % - apparent_range:
    %   - (1, 1) double.
    %   - The apparent range to the satellite.
    %   - Unit meters.
    % - apparent_range_rate:
    %   - (1, 1) double.
    %   - The apparent range rate to the satellite.
    %   - Unit meters / second.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th edition
    arguments(Input)
        satellite_position_at_rx_gcrf(3, 1) double {mustBeReal, mustBeFinite}
        satellite_velocity_at_rx_gcrf(3, 1) double {mustBeReal, mustBeFinite}
        station_position_at_rx_gcrf(3, 1) double {mustBeReal, mustBeFinite}
        station_velocity_at_rx_gcrf(3, 1) double {mustBeReal, mustBeFinite}
        options.threshold(1, 1) double {mustBeReal, mustBeFinite, mustBePositive} = 0.1 / Constants.SPEED_OF_LIGHT;
    end

    arguments(Output)
        apparent_range(1, 1) double {mustBeReal, mustBeFinite}
        apparent_range_rate(1, 1) double {mustBeReal, mustBeFinite}
        light_time(1, 1) double {mustBeReal, mustBeFinite}
    end

    % Assume constant velocity
    satellite_velocity_at_tx_gcrf = satellite_velocity_at_rx_gcrf; 

    % Approximate the initial light time
    satellite_position_at_tx_gcrf = satellite_position_at_rx_gcrf;
    light_time                    = vecnorm(station_position_at_rx_gcrf - satellite_position_at_tx_gcrf, 2, 1) / Constants.SPEED_OF_LIGHT;

    % Now iterate to find the satellite's position at transmission.
    threshold = options.threshold;
    while true
        satellite_position_at_tx_gcrf = satellite_position_at_rx_gcrf - satellite_velocity_at_rx_gcrf * light_time;

        new_range = vecnorm(station_position_at_rx_gcrf - satellite_position_at_tx_gcrf, 2, 1);
        new_light_time = new_range / Constants.SPEED_OF_LIGHT;

        if abs(new_light_time - light_time) < threshold
            break;
        else
            light_time = new_light_time;
        end
    end

    % Now evaluate the apparent range and range-rate
    apparent_range_vector = station_position_at_rx_gcrf - satellite_position_at_tx_gcrf;
    apparent_range        = vecnorm(apparent_range_vector, 2, 1);

    apparent_range_rate = dot(...
        (station_velocity_at_rx_gcrf - satellite_velocity_at_tx_gcrf), ...
        (apparent_range_vector ./ apparent_range));
end
