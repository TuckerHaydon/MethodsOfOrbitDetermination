%% Preamble
clc; clear all; close all;

%% Given parameters
gravitational_parameter = 3.986004415e14 * (Units.meters)^2 / (Units.seconds)^2;

%% Load the data
file_name = "1yr.180s.ECI.txt";
T = readtable(file_name, ...
    "FileType", "delimitedtext", ...
    "Delimiter", " ", ...
    "LeadingDelimitersRule", "ignore", ...
    "ConsecutiveDelimitersRule", "join", ...
    "NumHeaderLines", 0);

time                 = table2array(T(:, 1)); % seconds
time_utc_mjd         = table2array(T(:, 2)); % Julian days
position_history_eci = transpose(table2array(T(:, 3:5))); % meters
velocity_history_eci = transpose(table2array(T(:, 6:8))); % meters per second

%% Convert the ephemerides to orbital elements
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
eccentric_anomaly_history           = zeros(num_time, 1);
mean_anomaly_history                = zeros(num_time, 1);

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
        position_history_eci(:, idx), ...
        velocity_history_eci(:, idx), ...
        gravitational_parameter);

    eccentric_anomaly = EccentricAnomalyFromTrueAnomaly(true_anomaly, eccentricity);
    mean_anomaly      = MeanAnomalyFromEccentricAnomaly(eccentric_anomaly, eccentricity);

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
    eccentric_anomaly_history(idx)              = eccentric_anomaly;
    mean_anomaly_history(idx)                   = mean_anomaly;
end

%% Calculate the revolution times
% A revolution is defined as the time between two north-bound crossings of the equator (time when the satellite is at
% the line of nodes). This is when the true_anomaly is equal to 2*pi - argument_of_periapsis. Alternatively, this is
% when the ECI z-position crosses from negative to positive.

% This returns this index of the sample immediately before crossing 0.
idxs_immediately_before_crossing = find(diff(sign(position_history_eci(3, :))) > 0);
idxs_immediately_before_crossing = reshape(idxs_immediately_before_crossing, [], 1);
idxs_immediately_after_crossing  = idxs_immediately_before_crossing + 1;

num_crossings  = numel(idxs_immediately_before_crossing);
approx_crossing_times = zeros(num_crossings, 1);

for idx = 1:num_crossings
    % These indices capture the two points before zero crossing and the two points after.
    search_idxs = (-1:1:2) + idxs_immediately_before_crossing(idx);
    % Interpolating with pchip which does a shape-preserving piecewise cubic with 4 neighboring points.
    % In the future, might set up a least squares problem with the position/velocity constraints of the two points
    % surrounding the crossing.
    approx_crossing_times(idx) = interp1( ...
        position_history_eci(3, search_idxs), ...
        time(search_idxs), ...
        0, ...
        "pchip");
end

% Compute the starting and ending indices of each revoluton
mean_revolution_time     = mean(diff(approx_crossing_times));
num_full_num_revolutions = num_crossings - 1;
num_partial_revolutions  = num_full_num_revolutions + 2;

revolution_start_idxs = [1; idxs_immediately_after_crossing];
revolution_end_idxs   = [idxs_immediately_before_crossing; numel(time)];

revolution_start_times = time(revolution_start_idxs);
revolution_end_times = time(revolution_end_idxs);

approx_midpoint_revolution_times = (time(revolution_end_idxs) + time(revolution_start_idxs))/2;

%% Now compute the mean orbital elements
mean_semimajor_axis_history              = zeros(num_partial_revolutions, 1);
mean_eccentricity_history                = zeros(num_partial_revolutions, 1);
mean_inclination_history                 = zeros(num_partial_revolutions, 1);
mean_longitude_of_ascending_node_history = zeros(num_partial_revolutions, 1);
mean_argument_of_periapsis_history       = zeros(num_partial_revolutions, 1);
mean_mean_anomaly_history                = zeros(num_partial_revolutions, 1);

for revolution_idx = 1:num_partial_revolutions
    mean_semimajor_axis_history(revolution_idx)              = mean(semimajor_axis_history(              revolution_start_idxs(revolution_idx):revolution_end_idxs(revolution_idx)));
    mean_eccentricity_history(revolution_idx)                = mean(eccentricity_history(                revolution_start_idxs(revolution_idx):revolution_end_idxs(revolution_idx)));
    mean_inclination_history(revolution_idx)                 = mean(inclination_history(                 revolution_start_idxs(revolution_idx):revolution_end_idxs(revolution_idx)));
    mean_longitude_of_ascending_node_history(revolution_idx) = mean(longitude_of_ascending_node_history( revolution_start_idxs(revolution_idx):revolution_end_idxs(revolution_idx)));
    mean_argument_of_periapsis_history(revolution_idx)       = mean(argument_of_periapsis_history(       revolution_start_idxs(revolution_idx):revolution_end_idxs(revolution_idx)));
    mean_mean_anomaly_history(revolution_idx)                = mean(mean_anomaly_history(                revolution_start_idxs(revolution_idx):revolution_end_idxs(revolution_idx)));
end

%% Plot the mean and osculating elements
% Semimajor axis
figure('Position', [100, 100, 1500, 750]);
hold on;
plot( ...
    approx_midpoint_revolution_times ./ Units.days, ...
    mean_semimajor_axis_history ./ Units.kilometers, ...
    "-", ...
    "LineWidth", 3.0, ...
    "Color", "#0072BD");
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Time [days]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("a [km]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("Mean Semimajor Axis", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Eccentricity
figure('Position', [100, 100, 1500, 750]);
hold on;
plot( ...
    approx_midpoint_revolution_times ./ Units.days, ...
    mean_eccentricity_history, ...
    "-", ...
    "LineWidth", 3.0, ...
    "Color", "#0072BD");
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Time [days]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("e [unitless]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("Mean Eccentricity", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Inclination
figure('Position', [100, 100, 1500, 750]);
hold on;
plot( ...
    approx_midpoint_revolution_times ./ Units.days, ...
    rad2deg(mean_inclination_history), ...
    "-", ...
    "LineWidth", 3.0, ...
    "Color", "#0072BD");
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Time [days]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("I [deg]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("Mean Inclination", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Argument of periapsis
figure('Position', [100, 100, 1500, 750]);
hold on;
plot( ...
    approx_midpoint_revolution_times ./ Units.days, ...
    rad2deg(mean_argument_of_periapsis_history), ...
    "-", ...
    "LineWidth", 3.0, ...
    "Color", "#0072BD");
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Time [days]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("\omega [deg]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("Mean Argument of Periapsis", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Longitude of the ascending node
figure('Position', [100, 100, 1500, 750]);
hold on;
h2 = plot( ...
    approx_midpoint_revolution_times ./ Units.days, ...
    rad2deg(mean_longitude_of_ascending_node_history), ...
    "-", ...
    "LineWidth", 3.0, ...
    "Color", "#0072BD");
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Time [days]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("\Omega [deg]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("Mean Longitude of the Ascending Node", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

% Mean Anomaly
figure('Position', [100, 100, 1500, 750]);
hold on;
plot( ...
    approx_midpoint_revolution_times ./ Units.days, ...
    rad2deg(mean_mean_anomaly_history), ...
    "-", ...
    "LineWidth", 3.0, ...
    "Color", "#0072BD");
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Time [days]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
ylabel("M [deg]", "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
title("Mean Mean Anomaly", "FontSize", 18, "FontWeight", "Bold", "FontName", "Times New Roman");
axis tight;
ylim([-0.1, +0.1] .* diff(ylim()) + ylim());

%% Polyfit analysis
p = polyfit( ...
    approx_midpoint_revolution_times ./ Units.days, ...
    mean_semimajor_axis_history ./ Units.kilometers, ...
    1);
fprintf("Semimajor Axis: %0.3f km / day \n", p(1));

p = polyfit( ...
    approx_midpoint_revolution_times ./ Units.days, ...
    rad2deg(mean_inclination_history), ...
    1);
fprintf("Inclination: %0.3e deg / day \n", p(1));

p = polyfit( ...
    approx_midpoint_revolution_times ./ Units.days, ...
    rad2deg(mean_longitude_of_ascending_node_history), ...
    1);
fprintf("Longitude of Ascending Node: %0.3f deg / day \n", p(1));

idx = find(approx_midpoint_revolution_times ./ Units.days < 50, 1, 'last');
p = polyfit( ...
    approx_midpoint_revolution_times(1:idx) ./ Units.days, ...
    rad2deg(mean_argument_of_periapsis_history(1:idx)), ...
    1);
fprintf("Argument of Periapsis: %0.3f deg / day \n", p(1));

idx = find(approx_midpoint_revolution_times ./ Units.days < 50, 1, 'last');
p = polyfit( ...
    approx_midpoint_revolution_times(1:idx) ./ Units.days, ...
    rad2deg(mean_mean_anomaly_history(1:idx)), ...
    1);
fprintf("Mean Anomaly: %0.3f deg / day \n", p(1));