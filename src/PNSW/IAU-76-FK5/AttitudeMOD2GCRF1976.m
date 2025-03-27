function [C_mod2gcrf] = AttitudeMOD2GCRF1976( ...
        zeta_deg, ...
        theta_deg, ...
        z_deg)
    % Computes a coordinate transformation matrix that transforms a vector expressed in mean of date (MOD) frame to a
    % vector expressed in the geocentric celestial reference frame (GCRF). This transform accounts for the precession of
    % the Earth's pole and uses the IAU 1976 model.
    % 
    % Requires:
    % - zeta_deg:
    %   - (:, 1) double.
    %   - R3 precession rotation angle.
    %   - Unit degrees.
    % - theta_deg:
    %   - (:, 1) double.
    %   - R1 precession rotation angle.
    %   - Unit degrees.
    % - z_deg:
    %   - (:, 1) double.
    %   - R3 precession rotation angle.
    %   - Unit degrees.
    % 
    % Returns:
    % - C_mod2gcrf:
    %   - (3, 3, :) double.
    %   - The corresponding coordinate transform matrices from PEF to TOD.
    %
    % References:
    % - Vallado, Fundamental of Astrodynamics and applications, 5th edition
    arguments(Input)
        zeta_deg(:, 1) double {mustBeReal, mustBeFinite}
        theta_deg(:, 1) double {mustBeReal, mustBeFinite}
        z_deg(:, 1) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        C_mod2gcrf(3, 3, :) double {mustBeReal, mustBeFinite}
    end

    assert(all(size(zeta_deg) == size(theta_deg)), ...
        "Inputs must be same size.");
    assert(all(size(zeta_deg) == size(z_deg)), ...
        "Inputs must be same size.");

    % Vallado, 3-92
    C_mod2gcrf = pagemtimes( ...
        pagemtimes(ValladoROT3(zeta_deg), ValladoROT2(-theta_deg)), ...
        ValladoROT3(z_deg));
end