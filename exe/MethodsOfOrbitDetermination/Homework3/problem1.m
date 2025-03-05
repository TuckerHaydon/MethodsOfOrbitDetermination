%% Preamble
clc; clear all; close all;

%% Helper function
function [d_dt_state] = DynamicsFunction( ...
        time, ...
        state)
    arguments(Input)
        time(1, 1) double {mustBeReal, mustBeFinite}
        state(4, 1) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        d_dt_state(4, 1) double {mustBeReal, mustBeFinite}
    end

    % Rename
    x  = state(1);
    y  = state(2);
    vx = state(3);
    vy = state(4);
    r  = norm([x, y]);
    r3 = r^3;

    % Dynamics
    x_dot = vx;
    y_dot = vy;
    vx_dot = -x / r3;
    vy_dot = -y / r3;

    % Reform
    d_dt_state = [
        x_dot;
        y_dot;
        vx_dot;
        vy_dot;
    ];
end

function [F] = LinearizedDynamicsFunction( ...
        time, ...
        state)
    arguments(Input)
        time(1, 1) double {mustBeReal, mustBeFinite}
        state(4, 1) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        F(4, 4) double {mustBeReal, mustBeFinite}
    end

    % Rename
    x  = state(1);
    y  = state(2);
    vx = state(3);
    vy = state(4);
    r  = norm([x, y]);
    r5 = r^5;

    % Linearized Dynamics
    F11 = (2 * x^2 - y^2) / r5;
    F12 = 3 * x * y / r5;
    F21 = 3 * x * y / r5;
    F22 = (2 * y^2 - x^2) / r5;

    F = [
      0,   0    1, 0;
      0,   0    0, 1;
      F11, F12, 0, 0;
      F21, F22, 0, 0;
    ];
end

%% Double check 
% x = transpose([1, 0, 0, 1]);
% dx = 1e-4;
% 
% f_px = DynamicsFunction(0, x + [+dx; 0; 0; 0]);
% f_mx = DynamicsFunction(0, x + [-dx; 0; 0; 0]);
% 
% f_py = DynamicsFunction(0, x + [0; +dx; 0; 0]);
% f_my = DynamicsFunction(0, x + [0; -dx; 0; 0]);
% 
% f_pvx = DynamicsFunction(0, x + [0; 0; +dx; 0]);
% f_mvx = DynamicsFunction(0, x + [0; 0; -dx; 0]);
% 
% f_pvy = DynamicsFunction(0, x + [0; 0; 0; +dx]);
% f_mvy = DynamicsFunction(0, x + [0; 0; 0; -dx]);
% 
% F_approx = [
%     (f_px  - f_mx) / (2 * dx), ...
%     (f_py  - f_my) / (2 * dx), ...
%     (f_pvx - f_mvx) / (2 * dx), ...
%     (f_pvy - f_mvy) / (2 * dx)
% ];
% 
% F_exact = LinearizedDynamicsFunction(0, x);


%% Now numerically integrate this system
solution_times = (0:1:10) * 10;
initial_time   = solution_times(1);
final_time     = solution_times(end);

initial_true_state = transpose([1, 0, 0, 1]);

num_time     = numel(solution_times);
num_states   = numel(initial_true_state);

ode45_options = odeset('reltol', 3e-14, 'abstol', 1e-16);
[true_state_history, true_error_state_transition_matrix_history] = ...
    NumericallySolveTimeVaryingSystem( ...
        solution_times, ...
        initial_true_state, ...
        @DynamicsFunction, ...
        @LinearizedDynamicsFunction, ...
        ode45_options);

% Check
expected_true_state_1 = transpose([-0.839071529, -0.544021111, 0.544021111, -0.839071529]);
expected_true_state_10 = transpose([0.862318872, -0.506365641, 0.506365641, 0.862318872]);

assert(all(abs(true_state_history(:, 1+1) - expected_true_state_1) < 1e-9));
assert(all(abs(true_state_history(:, 10+1) - expected_true_state_10) < 1e-9));

% Report
fprintf("X[t_{1} ] = %0.3f\n", true_state_history(:, 1+1));
fprintf("X[t_{10}] = %0.3f\n", true_state_history(:, 1+10));

%% Now perturb the initial state
% true = nominal + error ---> nominal = true - error
initial_error_state   = transpose([1e-6, -1e-6, 1e-6, 1e-6]);
initial_nominal_state = initial_true_state - initial_error_state;

[nominal_state_history, nominal_error_state_transition_matrix_history] = ...
    NumericallySolveTimeVaryingSystem( ...
        solution_times, ...
        initial_nominal_state, ...
        @DynamicsFunction, ...
        @LinearizedDynamicsFunction, ...
        ode45_options);

solution_error = true_state_history - nominal_state_history;

true_final_solution_error    = true_error_state_transition_matrix_history(:, :, end)    * initial_error_state;
nominal_final_solution_error = nominal_error_state_transition_matrix_history(:, :, end) * initial_error_state;

% Check
expected_nominal_state_1 = transpose([-0.839031098, -0.544071486, 0.544076120, -0.839041244]);
expected_nominal_state_10 = transpose([0.862623360, -0.505843963, 0.505845689, 0.862623303]);

expected_state_transition_matrix_1 = [
    -19.2963174705, -1.0005919528, -1.5446240948, -20.5922746780;
    +24.5395368984, +2.5430400375, +3.3820224390, +24.9959638293;
    -26.6284485803, -1.2470410802, -2.0860289935, -27.5413748340;
    -15.0754226454, -1.4570972848, -2.0011442064, -14.6674122500;  
];

expected_state_transition_matrix_10 = [
    -151.2840323254, -0.0696433460, -0.5751839913, -152.5394552874;
    -260.2345144322, +0.8812356066, +0.0191322895, -260.6700884451;
    +259.1544475393, +0.3746434528, +1.2367484371, +260.0263802508;
    -152.1279107642, +0.3667128574, -0.1388295703, -151.6392131624;
];

assert(all(abs(nominal_state_history(:, 1+1)  - expected_nominal_state_1) < 1e-9));
assert(all(abs(nominal_state_history(:, 10+1) - expected_nominal_state_10) < 1e-9));

assert(all(abs(nominal_error_state_transition_matrix_history(:, :, 1+1) - expected_state_transition_matrix_1) < 1e-9, "all"));
assert(all(abs(nominal_error_state_transition_matrix_history(:, :, 10+1) - expected_state_transition_matrix_10) < 1e-9, "all"));


%% Symplectic matrix
J = [
    zeros(2, 2), eye(2, 2);
    -eye(2, 2), zeros(2, 2);
];

phi = nominal_error_state_transition_matrix_history(:, :, end);
phi_inv = inv(phi);
phi_inv_alt = -transpose(J * phi * J);
phi_inv_diff = phi_inv_alt - phi_inv;

%% Batch least squares
y = transpose([1, 2, 1]);
inv_R = diag([2, 1, 1]);
H = transpose([1, 1, 1]);

x_bar = 2;
inv_P_bar = 2;

cost_function = @(x) ...
    0.5 * transpose(y - H * x) * inv_R * (y - H * x) + ...
    0.5 * transpose(x_bar - x) * inv_P_bar * (x_bar - x);

x_hat = (transpose(H) * inv_R * H + inv_P_bar) \  (transpose(H) * inv_R * y + inv_P_bar * x_bar);

% Check that x_hat minimizes cost
cost_function(x_hat)
cost_function(x_hat + 0.01)
cost_function(x_hat - 0.01)

epsilon_hat = y - H * x_hat;