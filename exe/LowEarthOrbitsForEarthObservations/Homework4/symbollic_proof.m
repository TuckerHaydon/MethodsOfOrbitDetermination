clc; clear all; close all;

I         = sym("I", ["real"]);
Omega     = sym("Omega", ["real"]);
epsilon   = sym("epsilon", ["real"]);
Omega_sun = sym("Omega_sun", ["real"]);

% Compute e_hat_h
Rx = SymbollicValladoROT1(I);
Rz = SymbollicValladoROT3(Omega);

e_hat_h = transpose(Rz) * transpose(Rx) * [0; 0; 1];

% Compute e_hat_sun
Rx_sun = SymbollicValladoROT1(epsilon);
Rz_sun = SymbollicValladoROT3(Omega_sun);

e_hat_sun = transpose(Rx_sun) * transpose(Rz_sun) * [1; 0; 0];

sin_beta_prime = dot(e_hat_h, e_hat_sun);
sin_beta_prime_alt = sin(epsilon) * cos(I) * sin(Omega_sun) - cos(epsilon) * sin(I) * cos(Omega) * sin(Omega_sun) + sin(I) * sin(Omega) * cos(Omega_sun);

function [C] = SymbollicValladoROT1(theta_deg)
    % Reshape to row vector.
    theta_deg = reshape(theta_deg, 1, []);

    c     = cos(theta_deg);
    s     = sin(theta_deg);
    zero_ = zeros(size(c));
    one_  = ones(size(c));

    % Vallao, 3-15
    % Unrolled version
    C = [
        one_;
        zero_;
        zero_;
        zero_;
        c;
       -s;
        zero_;
        s;
        c;
    ];

    % Reshape to 3D array.
    C = reshape(C, 3, 3, []);
end

function [C] = SymbollicValladoROT3(theta_deg)
    % Reshape to row vector.
    theta_deg = reshape(theta_deg, 1, []);
    
    c     = cos(theta_deg);
    s     = sin(theta_deg);
    zero_ = zeros(size(c));
    one_  = ones(size(c));

    % Vallao, 3-15
    % Unrolled version
    C = [
        c;
       -s;
        zero_;
        s;
        c;
        zero_;
        zero_;
        zero_;
        one_;
    ];

    % Reshape to 3D array.
    C = reshape(C, 3, 3, []);
end