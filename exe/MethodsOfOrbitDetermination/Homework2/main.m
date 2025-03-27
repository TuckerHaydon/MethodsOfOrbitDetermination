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

earth_gravitational_parameter  = 398600.4 * (Units.kilometers)^3 / (Units.seconds)^2;
earth_dynamical_form_factor    = 0.00108248; % Unitless
earth_semimajor_axis           = 6378.145 * Units.kilometers;

%% Question 1: Numerically integrate the equations of motion using the two-body J2 gravitation model
% Let x = [position_eci, velocity_eci] be the state of the system.
% The system dynamics are:
%   - d/dt position_eci = velocity_eci;
%   - d/dt velocity_eci = gravitation_eci
dynamics_fun = @(t, x) [
    x(4:6); 
    ECIJ2Gravitation( ...
        x(1:3), ...
        earth_gravitational_parameter, ...
        earth_dynamical_form_factor, ...
        earth_semimajor_axis);
];

% Simulate for 1 day with 20-second increments
dt = 20 * Units.seconds;

time = colon(0, dt, 1 * Units.days);
time = reshape(time, [], 1);

% Initial state
x0 = [
    position_eci; 
    velocity_eci;
];

% RK4 integrator
ode45_options = odeset('reltol', 3e-14, 'abstol', 1e-16);
[~, x] = ode45(dynamics_fun, time, x0, ode45_options);

eci_position_history = transpose(x(:, 1:3));
eci_velocity_history = transpose(x(:, 4:6));

%% Check against provided values
final_position = transpose(x(end, :));
expected_final_position = [
    -5751.49900721589 * Units.kilometers;
    4721.14371040552  * Units.kilometers;
    2046.03583664311  * Units.kilometers;
    -0.797658631074   * Units.kilometers / Units.seconds;
    -3.656513108387   * Units.kilometers / Units.seconds;
    6.139612016678    * Units.kilometers / Units.seconds;
];

absolute_error = abs(final_position - expected_final_position);
relative_error = abs((final_position - expected_final_position) ./ expected_final_position);

assert(all(absolute_error(1:3) < 1e-3, 'all')); % mm accurate
assert(all(absolute_error(4:6) < 1e-6, 'all')); % um/s accuract
assert(all(relative_error < 1e-9, 'all'));

%% Convert the state vector to orbital elements
num_time = numel(time);

semiparameter_history               = zeros(num_time, 1);
semimajor_axis_history              = zeros(num_time, 1);
eccentricity_history                = zeros(num_time, 1);
inclination_history                 = zeros(num_time, 1);
longitude_of_ascending_node_history = zeros(num_time, 1);
argument_of_periapsis_history       = zeros(num_time, 1);
true_anomaly_history                = zeros(num_time, 1);
argument_of_latitude_history        = zeros(num_time, 1);
true_longitude_history              = zeros(num_time, 1);
true_longitude_of_periapsis_history = zeros(num_time, 1);

for idx = 1:num_time
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
        eci_position_history(:, idx), ...
        eci_velocity_history(:, idx), ...
        earth_gravitational_parameter);

    semiparameter_history(idx)                  = semiparameter;
    semimajor_axis_history(idx)                 = semimajor_axis;
    eccentricity_history(idx)                   = eccentricity;
    inclination_history(idx)                    = inclination;
    longitude_of_ascending_node_history(idx)    = longitude_of_ascending_node;
    argument_of_periapsis_history(idx)          = argument_of_periapsis;
    true_anomaly_history(idx)                   = true_anomaly;
    argument_of_latitude_history(idx)           = argument_of_latitude;
    true_longitude_history(idx)                 = true_longitude;
    true_longitude_of_periapsis_history(idx)    = true_longitude_of_periapsis;
end

eccentric_anomaly_history = EccentricAnomalyFromTrueAnomaly(true_anomaly_history, eccentricity_history);
mean_anomaly_history      = MeanAnomalyFromEccentricAnomaly(eccentric_anomaly_history, eccentricity_history);
orbital_period_history    = OrbitalPeriodFromSemimajorAxis(semimajor_axis_history, eccentricity_history, earth_gravitational_parameter);

%% Calculate the time since perigee
time_since_periapsis_history = (mean_anomaly_history / (2*pi)) .* orbital_period_history;

%% Plotting
% Semimajor axis
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, semimajor_axis_history, "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Semimajor Axis [m]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Semimajor Axis vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Eccentricity
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, eccentricity_history, "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Eccentricity [unitless]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Eccentricity vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Inclination
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, rad2deg(inclination_history), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Inclination [deg]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Inclination vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Longitude of Ascending Node
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, rad2deg(longitude_of_ascending_node_history), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Longitude of Ascending Node [deg]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Longitude of Ascending Node vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Argument of Periapsis
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, rad2deg(argument_of_periapsis_history), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Argument of Periapsis [deg]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Argument of Periapsis vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Time since periapsis
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, time_since_periapsis_history / Units.minutes, "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Time Since Periapsis [min]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Mean Time Since Periapsis vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

%% Compute expected orbital period
initial_orbital_period = OrbitalPeriodFromSemimajorAxis( ...
    semimajor_axis_history(1), ...
    eccentricity_history(1), ...
    earth_gravitational_parameter);

unwrapped_true_anomaly_history  = unwrap(true_anomaly_history);
num_completed_orbits            = floor(unwrapped_true_anomaly_history(end) / (2*pi));
orbit_completion_times          = interp1(unwrapped_true_anomaly_history, time, 2*pi*(1:num_completed_orbits));
orbit_periods                   = diff([0, orbit_completion_times]);

figure("Position", [100, 100, 1000, 750]);
hold on;
h1 = plot(1:num_completed_orbits, orbit_periods, "LineWidth", 3.0);
h2 = plot(1:num_completed_orbits, initial_orbital_period * ones(size(1:num_completed_orbits)), "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Completed Orbit Number", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Orbital Period [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Orbital Period Over Times", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
legend([h1, h2], "Actual Orbital Period", "Initial Expected Orbital Period", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

%% Frequency analysis of orbital elements
% time_per_orbit = time / orbit_periods(1);
% fs = 1.0/(time_per_orbit(2) - time_per_orbit(1));                               
% t = time_per_orbit;   
% 
% % x = argument_of_periapsis_history;
% % y = fft(argument_of_periapsis_history);
% 
% x = argument_of_periapsis_history - mean(argument_of_periapsis_history);
% y = fft(x);
% 
% % x = cos(4 * time_per_orbit) + cos(2 * time_per_orbit) + 5;
% % y = fft(x);
% 
% n = length(x);          % number of samples
% f = (0:n-1)*(fs/n);     % frequency range
% power = abs(y).^2/n;    % power of the DFT
% 
% y0 = fftshift(y);         % shift y values
% f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
% power0 = abs(y0).^2/n;    % 0-centered power
% 
% figure();
% hold on;
% plot(f0,power0, 'k-')
% xlabel('Frequency [Hz]')
% ylabel('Power')
% yscale('log')
% 
% [pks, locs] = findpeaks(power0);
% [~, I] = sort(pks);
% pks = pks(I);
% locs = locs(I);
% 
% % [~,I] = maxk(pks, 5);
% % pks = pks(I);
% % locs = locs(I);
% 
% freqs = f0(locs);
% % freqs(freqs < 0) = [];
% periods = 1./freqs;
% per_orbit = periods ./ 6740;
% 
% figure();
% stem(per_orbit, pks);

%% Compute the specific energy
U = ECIJ2GravitationalPotential( ...
        eci_position_history, ...
        earth_gravitational_parameter, ...
        earth_dynamical_form_factor, ...
        earth_semimajor_axis);

E = reshape((vecnorm(eci_velocity_history, 2, 1).^2) / 2, [], 1) - U;
dE = E - E(1);

figure("Position", [100, 100, 1000, 750]);
hold on;
plot(time ./ Units.hours, dE, "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Specific Energy [ J / kg ]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Specific Energy Relative to Epoch vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

%% Compute the specific angular momentum
h_vec_history = SpecificAngularMomentumFromStateVector( ...
    eci_position_history, ...
    eci_velocity_history);

dh_k = reshape(h_vec_history(3, :) - h_vec_history(3, 1), [], 1);

figure("Position", [100, 100, 1000, 750]);
hold on;
plot(time ./ Units.hours, dh_k, "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Specific Angular Momentum [ m^2 / s ]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Specific Z-Axis Angular Momentum Relative to Epoch vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

%% Integrate with drag
coefficient_of_drag     = 2.0; 
effective_area          = 3.6 * (Units.meters)^2;
mass                    = 1350 * Units.kilograms;
air_density             = 4e-13 * Units.kilograms / (Units.meters)^3;
reference_altitude      = 7298.145 * Units.kilometers;
altitude_rate           = 200.0 * Units.kilometers;
earth_rotation_rate     = 7.29211585530066e-5 * Units.radians / Units.seconds;

% Let x = [position_eci, velocity_eci] be the state of the system.
% The system dynamics are:
%   - d/dt position_eci = velocity_eci;
%   - d/dt velocity_eci = gravitation_eci + drag_eci
dynamics_fun_with_drag = @(t, x) [
    x(4:6); 
    ECIJ2Gravitation( ...
        x(1:3), ...
        earth_gravitational_parameter, ...
        earth_dynamical_form_factor, ...
        earth_semimajor_axis) + ...
    DragFunction( ...
        x(1:3), ...
        x(4:6), ...
        coefficient_of_drag, ...
        effective_area, ...
        mass, ...
        air_density, ...
        reference_altitude, ...
        altitude_rate, ...
        earth_rotation_rate);
];

% Simulate for 1 day with 20-second increments
dt = 20 * Units.seconds;

time = colon(0, dt, 1 * Units.days);
time = reshape(time, [], 1);

% Initial state
x0 = [
    position_eci; 
    velocity_eci;
];

% RK4 integrator
ode45_options = odeset('reltol', 3e-14, 'abstol', 1e-16);
[~, x_with_drag] = ode45(dynamics_fun_with_drag, time, x0, ode45_options);

eci_position_history_with_drag = transpose(x_with_drag(:, 1:3));
eci_velocity_history_with_drag = transpose(x_with_drag(:, 4:6));

%% Check against provided answer
final_position_with_drag = transpose(x_with_drag(end, :));
expected_final_position_with_drag = [
    -5751.50585441435 * Units.kilometers;
    4721.10759954747  * Units.kilometers;
    2046.09502951715  * Units.kilometers;
    -0.797610370476   * Units.kilometers / Units.seconds;
    -3.656553079577   * Units.kilometers / Units.seconds;
    6.139595227675    * Units.kilometers / Units.seconds;
];

absolute_error_with_drag = abs(final_position_with_drag - expected_final_position_with_drag);
relative_error_with_drag = abs((final_position_with_drag - expected_final_position_with_drag) ./ expected_final_position_with_drag);

assert(all(absolute_error_with_drag(1:3) < 1e-3, 'all')); % mm accurate
assert(all(absolute_error_with_drag(4:6) < 1e-6, 'all')); % um/s accuract
assert(all(relative_error_with_drag < 1e-9, 'all'));

%% Convert the state vector to orbital elements
num_time = numel(time);

semiparameter_history_with_drag               = zeros(num_time, 1);
semimajor_axis_history_with_drag              = zeros(num_time, 1);
eccentricity_history_with_drag                = zeros(num_time, 1);
inclination_history_with_drag                 = zeros(num_time, 1);
longitude_of_ascending_node_history_with_drag = zeros(num_time, 1);
argument_of_periapsis_history_with_drag       = zeros(num_time, 1);
true_anomaly_history_with_drag                = zeros(num_time, 1);
argument_of_latitude_history_with_drag        = zeros(num_time, 1);
true_longitude_history_with_drag              = zeros(num_time, 1);
true_longitude_of_periapsis_history_with_drag = zeros(num_time, 1);

for idx = 1:num_time
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
        eci_position_history_with_drag(:, idx), ...
        eci_velocity_history_with_drag(:, idx), ...
        earth_gravitational_parameter);

    semiparameter_history_with_drag(idx)                  = semiparameter;
    semimajor_axis_history_with_drag(idx)                 = semimajor_axis;
    eccentricity_history_with_drag(idx)                   = eccentricity;
    inclination_history_with_drag(idx)                    = inclination;
    longitude_of_ascending_node_history_with_drag(idx)    = longitude_of_ascending_node;
    argument_of_periapsis_history_with_drag(idx)          = argument_of_periapsis;
    true_anomaly_history_with_drag(idx)                   = true_anomaly;
    argument_of_latitude_history_with_drag(idx)           = argument_of_latitude;
    true_longitude_history_with_drag(idx)                 = true_longitude;
    true_longitude_of_periapsis_history_with_drag(idx)    = true_longitude_of_periapsis;
end

eccentric_anomaly_history_with_drag = EccentricAnomalyFromTrueAnomaly(true_anomaly_history_with_drag, eccentricity_history_with_drag);
mean_anomaly_history_with_drag      = MeanAnomalyFromEccentricAnomaly(eccentric_anomaly_history_with_drag, eccentricity_history_with_drag);
orbital_period_history_with_drag    = OrbitalPeriodFromSemimajorAxis(semimajor_axis_history_with_drag, eccentricity_history_with_drag, earth_gravitational_parameter);

%% Calculate the time since perigee
time_since_periapsis_history_with_drag = (mean_anomaly_history_with_drag / (2*pi)) .* orbital_period_history_with_drag;

%% Plotting
% Semimajor axis
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, semimajor_axis_history_with_drag - semimajor_axis_history, "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Semimajor Axis Difference [m]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Semimajor Axis Difference with Drag vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Eccentricity
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, eccentricity_history_with_drag - eccentricity_history, "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Eccentricity Difference [unitless]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Eccentricity Difference with Drag vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Inclination
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, rad2deg(inclination_history_with_drag - inclination_history), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Inclination Difference [deg]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Inclination Difference with Drag vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Longitude of Ascending Node
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, rad2deg(longitude_of_ascending_node_history_with_drag - longitude_of_ascending_node_history), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Longitude of Ascending Node Difference [deg]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Longitude of Ascending Node Difference with Drag vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Argument of Periapsis
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, rad2deg(argument_of_periapsis_history_with_drag - argument_of_periapsis_history), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Argument of Periapsis Difference [deg]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Argument of Periapsis Difference with Drag vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Time since periapsis
figure("Position", [100, 100, 1000, 750]);
plot(time ./ Units.hours, (time_since_periapsis_history_with_drag - time_since_periapsis_history), "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Time Since Periapsis Difference [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Mean Time Since Periapsis Difference vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());


%% Compute expected orbital period
initial_orbital_period_with_drag = OrbitalPeriodFromSemimajorAxis( ...
    semimajor_axis_history_with_drag(1), ...
    eccentricity_history_with_drag(1), ...
    earth_gravitational_parameter);

unwrapped_true_anomaly_history_with_drag    = unwrap(true_anomaly_history_with_drag);
num_completed_orbits_with_drag              = floor(unwrapped_true_anomaly_history_with_drag(end) / (2*pi));
orbit_completion_times_with_drag            = interp1(unwrapped_true_anomaly_history_with_drag, time, 2*pi*(1:num_completed_orbits));
orbit_periods_with_drag                     = diff([0, orbit_completion_times_with_drag]);

figure("Position", [100, 100, 1000, 750]);
hold on;
plot( ...
    1:num_completed_orbits_with_drag, ...
    orbit_periods_with_drag - orbit_periods, ...
    "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Completed Orbit Number", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Orbital Period Difference [s]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Orbital Period Difference with Drag Over Times", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

%% Compute the specific energy
U_with_drag = ECIJ2GravitationalPotential( ...
        eci_position_history_with_drag, ...
        earth_gravitational_parameter, ...
        earth_dynamical_form_factor, ...
        earth_semimajor_axis);

E_with_drag = reshape((vecnorm(eci_velocity_history_with_drag, 2, 1).^2) / 2, [], 1) - U_with_drag;
dE_with_drag = E_with_drag - E_with_drag(1);

figure("Position", [100, 100, 1000, 750]);
hold on;
h1 = plot(time ./ Units.hours, dE_with_drag, "k-", "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
xlabel("Time [hr]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
ylabel("Specific Energy [ J / kg ]", "FontSize", 14, "FontWeight", "Bold", "FontName", "TimesNewRoman");
title("Specific Energy Relative to Epoch with Drag vs Time", "FontSize", 18, "FontWeight", "Bold", "FontName", "TimesNewRoman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

%% Helper functions
function [acceleration_eci] = DragFunction( ...
        position_eci, ...
        velocity_eci, ...
        coefficient_of_drag, ...
        effective_area, ...
        mass, ...
        air_density, ...
        reference_altitude, ...
        altitude_rate, ...
        earth_rotation_rate)

    % Rename
    r_vec = position_eci;
    v_vec = velocity_eci;
    C_d = coefficient_of_drag;
    A = effective_area;
    m = mass;
    rho_0 = air_density;
    r_0 = reference_altitude;
    H = altitude_rate;
    theta_dot = earth_rotation_rate;
    
    r_mag = vecnorm(r_vec, 2, 1);

    rho_A = rho_0 .* exp(-(r_mag - r_0) ./ H);

    V_A_vec = v_vec + [
        theta_dot .* r_vec(2, :); 
        -theta_dot .* r_vec(1, :); 
        zeros(size(r_vec(1, :)));
    ];
    V_A = vecnorm(V_A_vec, 2, 1);

    % Compute the drag 
    acceleration_eci = -1/2 .* C_d .* A ./ m .* rho_A .* V_A .* V_A_vec;
end