%% Preamble
clc; clear all; close all;

%% Setup
mu = sym("mu", ["real", "positive"]);
J2 = sym("J2", ["real", "positive"]);
Re = sym("Re", ["real", "positive"]);
x = sym("x", [3, 1], "real");

r = norm(x);
sin_phi = x(3) / r;

U = (mu / r) * (1 - J2 * (Re/r)^2 * ((3/2) * sin_phi^2 - 1/2));
U_fun = matlabFunction(U, 'Vars', transpose([mu; J2; Re; x]));

%% Compute the Jacobian
jacobian_U = jacobian(U, x);

%% Compare against hand-calculated jacobian with Monte Carlo
u_hat_z = [0; 0; 1];
x_hat = x ./ r;
% hand_calculated_jacobian_U = -(mu / r^3) * (transpose(x) + (3 * J2 * Re^2 / r) * ((1/2 - 5/2 * sin_phi^2) * transpose(x_hat) + sin_phi * transpose(u_hat_z)));
hand_calculated_jacobian_U = -(mu / r^3) * (transpose(x) + (3/2) * J2 * (Re / r)^2 * [
    (1 - 5*(x(3)/r)^2) * x(1), ...
    (1 - 5*(x(3)/r)^2) * x(2), ...
    (3 - 5*(x(3)/r)^2) * x(3), ...
]);

jacobian_U_fun = matlabFunction(jacobian_U, 'Vars', transpose([mu; J2; Re; x]));
hand_calculated_jacobian_U_fun = matlabFunction(hand_calculated_jacobian_U, 'Vars', transpose([mu; J2; Re; x]));

mu_ = 398600.4*10^9;
J2_ = 0.00108248;
Re_ = 6378.145e3;

num_iter = 1e5;

differences = zeros(num_iter, 3);

for idx = 1:num_iter

    % https://math.libretexts.org/Courses/Mount_Royal_University/MATH_2200%3A_Calculus_for_Scientists_II/7%3A_Vector_Spaces/5.7%3A_Cylindrical_and_Spherical_Coordinates
    theta = 2 * pi * rand();
    phi = pi * rand();
    r = Re_;

    x_ = [
        r * sin(phi) * cos(theta);
        r * sin(phi) * sin(theta);
        r * cos(phi);
    ];

    analytical_g = jacobian_U_fun(mu_, J2_, Re_, x_(1), x_(2), x_(3));
    hand_calculated_g = hand_calculated_jacobian_U_fun(mu_, J2_, Re_, x_(1), x_(2), x_(3));

    differences(idx, :) = hand_calculated_g - analytical_g;
end

max(vecnorm(differences, 2, 2))