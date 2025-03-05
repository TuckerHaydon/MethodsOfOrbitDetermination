function [zeta_deg, theta_deg, z_deg] = ComputePrecessionParameters1976( ...
        julian_centuries_since_epoch)
    % Computes the Earth precession parameters using the IAU 76 model.
    % 
    % Requires:
    % - julian_centuries_since_epoch:
    %   - (:, :) double.
    %   - The number of julian centuries since the J2000.0 epoch.
    % 
    % Returns:
    % - zeta_deg:
    %   - (:, :) double.
    %   - R3 precession rotation angle.
    %   - Unit degrees.
    % - theta_deg:
    %   - (:, :) double.
    %   - R1 precession rotation angle.
    %   - Unit degrees.
    % - z_deg:
    %   - (:, :) double.
    %   - R3 precession rotation angle.
    %   - Unit degrees.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th edition
    arguments(Input)
        julian_centuries_since_epoch(:, :) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        zeta_deg(:, :) double {mustBeReal, mustBeFinite}
        theta_deg(:, :) double {mustBeReal, mustBeFinite}
        z_deg(:, :) double {mustBeReal, mustBeFinite}
    end

    % Rename
    TTT = julian_centuries_since_epoch;

    % Vallado, 3-91
    zeta_arcseconds  = 2306.2181 .* TTT + 0.30188 .* TTT.^2 + 0.017998 .* TTT.^3;
    theta_arcseconds = 2004.3109 .* TTT - 0.42665 .* TTT.^2 - 0.041833 .* TTT.^3;
    z_arcsecond      = 2306.2181 .* TTT + 1.09468 .* TTT.^2 + 0.018203 .* TTT.^3;

    zeta_deg  = zeta_arcseconds  * Conversions.ARCSECONDS_TO_DEGREES;
    theta_deg = theta_arcseconds * Conversions.ARCSECONDS_TO_DEGREES;
    z_deg     = z_arcsecond      * Conversions.ARCSECONDS_TO_DEGREES;
end