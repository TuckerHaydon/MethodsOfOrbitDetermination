%% Preamble
clc; clear all; close all;

%% Given
position_eci = [
    -2436.45;
    -2436.45;
    6891.037;
] * Units.kilometers;

velocity_eci = [
    +5.088611;
    -5.088611;
    +0.0;
] * Units.kilometers / Units.seconds;

gravitational_parameter = 398600.5 * (Units.kilometers)^3 / (Units.seconds)^2;

%% Question 1: Compute the orbital elements
[...
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
        gravitational_parameter);

%% Question 2: Compute the ECI vector
recalculated_semiparameter = SemiparameterFromSemimajorAxis( ...
    semimajor_axis, eccentricity);

[...
        recalculated_position_eci, ...
        recalculated_velocity_eci] = ...
    StateVectorFromOrbitalElements( ...
        recalculated_semiparameter, ...
        eccentricity, ...
        inclination, ...
        longitude_of_ascending_node, ...
        argument_of_periapsis, ...
        true_anomaly, ...
        argument_of_latitude, ...
        true_longitude, ...
        true_longitude_of_periapsis, ...
        gravitational_parameter);

%% Question 4: Numerically integrate the equations of motion using the two-body gravitation model
% Let x = [position_eci, velocity_eci] be the state of the system.
% The system dynamics are:
%   - d/dt position_eci = velocity_eci;
%   - d/dt velocity_eci = gravitation_eci
dynamics_fun = @(t, x) [x(4:6); ECISimpleTwoBodyGravitation(x(1:3), gravitational_parameter)];

% Want to simulate two full orbits.
orbital_period = OrbitalPeriodFromSemimajorAxis( ...
    semimajor_axis, ...
    eccentricity, ...
    gravitational_parameter);

dt = 20 * Units.seconds;

time = colon(0, dt, 2 * orbital_period);
time = reshape(time, [], 1);

% Initial state
x0 = [
    position_eci; 
    velocity_eci;
];

% RK4 integrator
ode45_options = odeset('reltol', 1e-6, 'abstol', 1e-6);
% ode45_options = odeset('RelTol', 1e-5, 'MaxStep', 1e-6);
% ode45_options = odeset();
[~, x] = ode45(dynamics_fun, time, x0, ode45_options);

position_eci_history = transpose(x(:, 1:3));
velocity_eci_history = transpose(x(:, 4:6));
acceleration_eci_history = ECISimpleTwoBodyGravitation( ...
    position_eci_history, ...
    gravitational_parameter);

specific_angular_momentum_history = SpecificAngularMomentumFromStateVector( ...
    position_eci_history, ...
    velocity_eci_history);

%% Question 4 Plotting
figure("Position", [100, 100, 900, 900]);
plot3( ...
    position_eci_history(1, :), ...
    position_eci_history(2, :), ...
    position_eci_history(3, :), ...
    "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("ECI X [m]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("ECI Y [m]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
zlabel("ECI Z [m]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Orbit", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

figure('Position', [100, 100, 900, 600]);
scatter3( ...
    specific_angular_momentum_history(1, :), ...
    specific_angular_momentum_history(2, :), ...
    specific_angular_momentum_history(3, :), ...
    "k.", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("X [m^2 / s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Y [m^2 / s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
zlabel("Z [m^2 / s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Specific Angular Momentum", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");


figure("Position",[100, 100, 1000, 1000]);
tlt = tiledlayout(3, 3);
title(tlt, "ECI Position, Velocity, & Acceleration", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(1);
plot(time, position_eci_history(1, :), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
% xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Position [m]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("X", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(4);
plot(time, velocity_eci_history(1, :), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
% xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Velocity [m/s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
% title("X", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(7);
plot(time, acceleration_eci_history(1, :), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Acceleration [m/s^2]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
% title("X", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(2);
plot(time, position_eci_history(2, :), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
% xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Position [m]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Y", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(5);
plot(time, velocity_eci_history(2, :), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
% xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Velocity [m/s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
% title("Y", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(8);
plot(time, acceleration_eci_history(2, :), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Acceleration [m/s^2]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
% title("Y", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(3);
plot(time, position_eci_history(3, :), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
% xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Position [m]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Z", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(6);
plot(time, velocity_eci_history(3, :), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
% xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Velocity [m/s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
% title("Z", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(9);
plot(time, acceleration_eci_history(3, :), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Acceleration [m/s^2]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
% title("Z", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");


%% Question 5: Compute the specific mechanical/potential energy as a function of time.
specific_kinetic_energy_history   = 0.5 * dot(velocity_eci_history, velocity_eci_history, 1);
specific_potential_energy_history = -gravitational_parameter ./ sqrt(dot(position_eci_history, position_eci_history, 1));

specific_mechanical_energy_history = SpecificMechanicalEnergyFromStateVector( ...
        position_eci_history, ...
        velocity_eci_history, ...
        gravitational_parameter);

% Plotting
figure("Position", [100, 100, 1500, 750]);
tlt = tiledlayout(2, 2);

nexttile(1);
plot(time, specific_kinetic_energy_history, "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Energy [m^2 / s^2]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Specific Kinetic Energy vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(3);
plot(time, specific_potential_energy_history, "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Energy [m^2 / s^2]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Specific Potential Energy vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(2);
plot(time, specific_mechanical_energy_history, "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Energy [m^2 / s^2]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Specific Mechanical Energy vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");

nexttile(4);
plot(time, specific_mechanical_energy_history - specific_mechanical_energy_history(1), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Energy Deviation [m^2 / s^2]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Specific Mechanical Energy Deviation vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
