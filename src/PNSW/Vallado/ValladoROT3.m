function [C] = ValladoROT3(theta_deg)
    % Elementary rotation about the z axis.
    %
    % Requires:
    % - theta_deg:
    %   - (:, 1) double.
    %   - Set of angles to rotate about the z axis.
    %   - Unit degrees.
    % 
    % Returns:
    % - C:
    %   - (3, 3, :) double.
    %   - Corresponding rotation matrices.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th edition
    arguments(Input)
        theta_deg(:, 1) {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        C(3, 3, :) {mustBeReal, mustBeFinite}
    end

    % Reshape to row vector.
    theta_deg = reshape(theta_deg, 1, []);
    
    c     = cosd(theta_deg);
    s     = sind(theta_deg);
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