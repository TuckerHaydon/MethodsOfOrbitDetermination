function [d_dt_state] = CubatureTransformContinuousTimeDynamicsFunction(...
        utc_time, ...
        state, ...
        state_parameters, ...
        time_parameters, ...
        gravity_parameters, ...
        drag_parameters, ...
        structure_parameters)
    arguments(Input)
        utc_time(1, 1) datetime
        state(:, :) double
        state_parameters(1, 1) struct
        time_parameters(1, 1) struct
        gravity_parameters(1, 1) struct
        drag_parameters(1, 1) struct
        structure_parameters(1, 1) struct
    end

    arguments(Output)
        d_dt_state(:, :) double
    end

    % State should be size [num_states, N]
    assert(size(state, 1) == state_parameters.num_states);

    % Unpack the state samples
    satellite_position_gcrf = state(1:3, :);
    satellite_velocity_gcrf = state(4:6, :);
    coefficient_of_drag     = state(7, :);

    % Push the state samples through the dynamics function
    [d_dt_satellite_position_gcrf, d_dt_satellite_velocity_gcrf, d_dt_coefficient_of_drag] = ...
        DynamicsFunction(...
            utc_time, ...
            satellite_position_gcrf, ...
            satellite_velocity_gcrf, ...
            coefficient_of_drag, ...
            time_parameters, ...
            gravity_parameters, ...
            drag_parameters, ...
            structure_parameters);

    % Pack the state samples
    d_dt_state = [
        d_dt_satellite_position_gcrf; 
        d_dt_satellite_velocity_gcrf; 
        reshape(d_dt_coefficient_of_drag, 1, []);
    ];

    % State should be size [num_states, N]
    assert(all(size(d_dt_state) == size(state), 'all'));
end