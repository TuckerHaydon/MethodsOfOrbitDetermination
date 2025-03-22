%% Preamble
clc; clear all; close all;

%% Problem 2
% Earth
% J2 = 0.001082;
% ae = 6378136.3;
% mu = 3.986004415e14;

% Mars
% https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
% https://en.wikipedia.org/wiki/Standard_gravitational_parameter
J2 = 1960.45e-6;
ae = 3396.2e3;
mu = 4.2828372e13;

inclinations = reshape(deg2rad(1:1:179), 1, []);

eccentricities = [
    0;
    0.030;
    0.0065;
    0.128;
    0.183;
    0.231;
    0;
];

periapsis_altitudes = [
    100;
    100;
    100;
    100;
    100;
    100;
    4000;
] .* 1e3;

apoapsis_altitudes = [
    100;
    500;
    1000;
    2000;
    3000;
    4000;
    4000;
] .* 1e3;


semimajor_axes = (apoapsis_altitudes + periapsis_altitudes + 2 .* ae) ./ 2;

mean_motions = sqrt(mu ./ semimajor_axes.^3);


% Rename
n_bar = mean_motions;
a = semimajor_axes;
apsidal_rate = -(3/4) .* n_bar .* (ae ./ a).^2 .* J2 .* (1 - 5 .* cos(inclinations).^2) ./ (1 - eccentricities.^2).^2;

figure("Position", [100, 100, 750, 500]);
h = plot(rad2deg(inclinations), rad2deg(apsidal_rate) * (60 * 60 * 24), "LineWidth", 3.0);
grid on;
grid minor;
set(gca, "GridAlpha", 0.75, "MinorGridAlpha", 0.5);
set(gca, "FontSize", 14, "FontWeight", "Bold", "FontName", "Times New Roman");
xlabel("Inclination [deg]", "FontSize", 14, "FontWeight", "bold", "FontName", "Times New Roman");
ylabel("Apsidal Rate [deg / day]", "FontSize", 14, "FontWeight", "bold", "FontName", "Times New Roman");
ylim([-10, 30]);
legend([h(1), h(2), h(3), h(4), h(5), h(6), h(7)], ...
    sprintf("[%0.3f, %d km, %d km]", eccentricities(1), periapsis_altitudes(1) / 1e3, apoapsis_altitudes(1) / 1e3), ...
    sprintf("[%0.3f, %d km, %d km]", eccentricities(2), periapsis_altitudes(2) / 1e3, apoapsis_altitudes(2) / 1e3), ...
    sprintf("[%0.3f, %d km, %d km]", eccentricities(3), periapsis_altitudes(3) / 1e3, apoapsis_altitudes(3) / 1e3), ...
    sprintf("[%0.3f, %d km, %d km]", eccentricities(4), periapsis_altitudes(4) / 1e3, apoapsis_altitudes(4) / 1e3), ...
    sprintf("[%0.3f, %d km, %d km]", eccentricities(5), periapsis_altitudes(5) / 1e3, apoapsis_altitudes(5) / 1e3), ...
    sprintf("[%0.3f, %d km, %d km]", eccentricities(6), periapsis_altitudes(6) / 1e3, apoapsis_altitudes(6) / 1e3), ...
    sprintf("[%0.3f, %d km, %d km]", eccentricities(7), periapsis_altitudes(7) / 1e3, apoapsis_altitudes(7) / 1e3), ...
"FontSize", 14, "FontWeight", "bold", "FontName", "Times New Roman", "Location", "North");
title("Apsidal Rate vs Inclination and Orbital Eccentricity", "FontSize", 18, "FontWeight", "bold", "FontName", "Times New Roman");

