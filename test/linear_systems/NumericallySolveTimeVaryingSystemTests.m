classdef NumericallySolveTimeVaryingSystemTests < matlab.unittest.TestCase
    methods (Test)
        function TestLinearTimeInvariantSystem(testCase)   

            A = [
                0,    1;
                -0.1, 0;
            ];

            initial_time = 0;
            initial_state = [
                0; 
                1;
            ];

            final_time = 10;
            expected_discrete_time_error_state_transition_matrix = expm(A * (final_time - initial_time));
            expected_final_state = expected_discrete_time_error_state_transition_matrix * initial_state;

            state_dynamics_function = @(t, x) A * x;
            linearized_state_dynamics_function = @(t, x) A;
            ode45_options = odeset('RelTol', 1e-12, "AbsTol", 1e-12);

            [final_state, discrete_time_error_state_transition_matrix] = ...
                NumericallySolveTimeVaryingSystem( ...
                    initial_time, ...
                    final_time, ...
                    initial_state, ...
                    state_dynamics_function, ...
                    linearized_state_dynamics_function, ...
                    ode45_options);

            testCase.verifyEqual( ...
                final_state, ...
                expected_final_state, ...
                'AbsTol', 1e-9); 
            testCase.verifyEqual( ...
                discrete_time_error_state_transition_matrix, ...
                expected_discrete_time_error_state_transition_matrix, ...
                'AbsTol', 1e-9); 
        end

    end
end