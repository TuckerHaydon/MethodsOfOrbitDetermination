clc; clear all; close all;

epsilon_bar = sym("epsilon_bar", ["real"]);
delta_psi   = sym("delta_psi", ["real"]);
epsilon     = sym("epsilon", ["real"]);

C_tod2mod = SymbollicValladoROT1(-epsilon_bar) * SymbollicValladoROT3(delta_psi) * SymbollicValladoROT1(epsilon);


function [C] = SymbollicValladoROT1(theta)
    % Reshape to row vector.
    theta = reshape(theta, 1, []);

    c     = cos(theta);
    s     = sin(theta);
    zero_ = zeros(size(c));
    one_  = ones(size(c));

    % Vallado, 3-15
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

function [C] = SymbollicValladoROT3(theta)
    % Reshape to row vector.
    theta = reshape(theta, 1, []);
    
    c     = cos(theta);
    s     = sin(theta);
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