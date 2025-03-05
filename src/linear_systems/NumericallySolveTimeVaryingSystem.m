function [state_history, error_state_transition_matrix_history] = ...
    NumericallySolveTimeVaryingSystem( ...
        solution_times, ...
        initial_state, ...
        state_dynamics_function, ...
        linearized_state_dynamics_function, ...
        ode45_options)
    % Numerically integrates a time-varying system and computes the final state and error state transition matrix.
    %
    % Suppose we are given a true time-varying system of the form
    %   d/dt x_star[t] = f[x_star[t]]
    % where x_star is the true state, t is time, and f[] is the dynamics function.
    % 
    % Define the relationship between the true, nominal, and error states as
    %   x_star[t] = x[t] + delta_x[t]
    % where x is the nominal state and delta_x is the true state. Consequently, 
    %   d/dt x_star[t] = d/dt x[t] + d/dt delta_x[t]
    % 
    % Expand the dynamics function to the first order about the nominal state x:
    %   f[x_star[t]] ~ f[x[t]] + F[x[t]] * (x_star[t] - x[t]) + HOT
    % where F[x[t]] is the partial derivative of the dynamics function w.r.t. the state evaluated at the nominal state, and
    % HOT are higher-order terms to be neglected.
    % 
    % Substituting into the state dynamics and simplifying produces
    %   d/dt x[t] + d/dt delta_x[t] = f[x[t]] + F[x[t]] * delta_x[t]
    % which, by the principle of superposition, can be written as the sum of the following two systems:
    %   d/dt x[t]       = f[x[t]]
    %   d/dt delta_x[t] = F[x[t]] * delta_x[t]
    %
    % The nominal system must be solved via numerical integration. The error system, too, could be solved with numerical
    % integration, however its solution could also be written as: 
    %   delta_x[t] =  phi[t, tau] delta_x[tau]
    % where phi[t, tau] is the discrete-time error state transition matrix. This matrix is required to propagate
    % the system's error covariance, so it is prudent to solve the error-state system by computing this matrix. This
    % matrix may be computed by numerically solving the following matrix differential equation: 
    %   d/dt phi[t, tau] = F[x[t]] phi[t, tau]
    % Note that this differential equation implicitly depends on the state x[t], hence why the state and error state
    % transition matrix must be numerically integrated together.
    % 
    % Requires:
    % Returns:
    % References:
    arguments(Input)
        solution_times(:, 1) double {mustBeReal, mustBeFinite}
        initial_state(:, 1) double {mustBeReal, mustBeFinite}
        state_dynamics_function(1, 1) function_handle
        linearized_state_dynamics_function(1, 1) function_handle
        ode45_options(1, 1) struct
    end

    arguments(Output)
        state_history(:, :) double {mustBeReal, mustBeFinite}
        error_state_transition_matrix_history(:, :, :) double {mustBeReal, mustBeFinite}
    end

    assert(all(diff(solution_times) > 0), "State must be propagated *forward* in time");

    num_solution_times = numel(solution_times);
    num_states         = numel(initial_state);

    % The initial state transition matrix is always identity.
    initial_state_transition_matrix = eye(num_states);

    % Form an augmented state. The augmented state is defined as the state stacked on top of the vectorized state
    % transition matrix.
    initial_augmented_state = [
        initial_state(:); 
        initial_state_transition_matrix(:);
    ];

    % Form the augmented dynamics function.
    augmented_dynamics_fun = @(time, augmented_state) AugmentedStateDynamicsFunction(...
        time, ...
        augmented_state, ...
        num_states, ...
        state_dynamics_function, ...
        linearized_state_dynamics_function);

    % Now numerically integrate with ode45.
    [~, augmented_state_solution] = ode45(...
        augmented_dynamics_fun, ...
        solution_times, ...
        initial_augmented_state, ...
        ode45_options);

    % Transpose the solution so its [num_states x num_times]
    augmented_state_solution = transpose(augmented_state_solution);

    % Extract the state and error state transition matrix from the augmented state.
    state_history                         = augmented_state_solution(1:num_states, :);
    error_state_transition_matrix_history = augmented_state_solution((num_states + 1):end, :);

    % Reshape the error state transition matrix history
    error_state_transition_matrix_history = reshape( ...
        error_state_transition_matrix_history, ...
        [num_states, num_states, num_solution_times]);

    % Output checks
    assert(all([num_states, num_solution_times] == size(state_history)));
    assert(all([num_states, num_states, num_solution_times] == size(error_state_transition_matrix_history)));

end

function [d_dt_augmented_state] = AugmentedStateDynamicsFunction(...
        time, ...
        augmented_state, ...
        num_states, ...
        state_dynamics_function, ...
        linearized_state_dynamics_function)

    % Expand the augmented state
    state = augmented_state(1:num_states);
    state_transition_matrix = reshape(augmented_state((num_states + 1):end), num_states, num_states);

    % d/dt x = f[x]
    d_dt_state = state_dynamics_function(time, state);

    % d/dt phi = F[x] * phi
    d_dt_state_transition_matrix = linearized_state_dynamics_function(time, state) * state_transition_matrix;

    % Now reform the augmented state
    d_dt_augmented_state = [
       d_dt_state;
       d_dt_state_transition_matrix(:);
    ];
end