function [C_itrf2pef] = AttitudeITRF2PEF1976(polar_motion_deg)
    % Computes the coordinate transform matrix that transforms a vector expressed in the internation terrestrial
    % reference frame (ITRF) to a vector expressed in the pseudo-Earth frame (PEF). This transform accounts for the
    % change in the instantaneous Earth rotation pole and is based on the IAU 1976 reduction.
    % 
    % Requires:
    % - polar_motion_deg:
    %   - (2, :) double.
    %   - The x and y polar motion angles.
    %   - Unit degrees.
    % 
    % Returns:
    % - C_itrf2pef:
    %   - (3, 3, :) double.
    %   - The coordinate transform matrix.
    % 
    % References:
    % - Vallado, Fundamental of Astrodynamics and applications
    % - Polar motion angles retrieved from: 
    %   - https://www.iers.org/IERS/EN/Publications/Bulletins/bulletins.html
    arguments(Input)
        polar_motion_deg(2, :) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        C_itrf2pef(3, 3, :) double {mustBeReal, mustBeFinite}
    end

    num_dims = size(polar_motion_deg, 2);

    % Rename
    x = polar_motion_deg(1, :);
    y = polar_motion_deg(2, :);

    % Vallado, Eq. 3-81
    % C_itrf2pef = [
    %     c[x],        0,          -s[x];
    %     s[x]s[y],    c[y],    c[x]s[y];
    %     s[x]c[y],   -s[x],    c[x]c[y];
    % ];

    % This is the vectorized version of above
    cx_   = cosd(x);
    sx_   = sind(x);
    cy_   = cosd(y);
    sy_   = sind(y);
    zero_ = zeros(1, num_dims);

    C_itrf2pef = [
        cx_;
        sx_ .* sy_;
        sx_ .* cy_;
        zero_;
        cy_;
        -sy_;
        -sx_;
        cx_ .* sy_;
        cx_ .* cy_;
    ];
    C_itrf2pef = reshape(C_itrf2pef, [3, 3, num_dims]);
end