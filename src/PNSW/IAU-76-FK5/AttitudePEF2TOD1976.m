function [C_pef2tod] = AttitudePEF2TOD1976(GAST_deg)
    % Computes a coordinate transformation matrix that transforms a vector expressed in the pseudo-Earth frame (PEF)
    % into a vector expressed in the true of date (TOD) frame. This transform accounts for the rotation of the Earth
    % about its pole and is based on the IAU-76 reduction.
    % 
    % Requires:
    % - GAST_deg:
    %   - (:, 1) double.
    %   - The Greenwich apparent sidereal time.
    %   - Unit degrees.
    % 
    % Returns:
    % - C_pef2tod:
    %   - (3, 3, :) double.
    %   - The corresponding coordinate transform matrices from PEF to TOD.
    %
    % References:
    % - Vallado, Fundamental of Astrodynamics and applications, 5th edition
    arguments(Input)
        GAST_deg(:, 1) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        C_pef2tod(3, 3, :) double {mustBeReal, mustBeFinite}
    end

    % Vallado, 3-83
    C_pef2tod = ValladoROT3(-GAST_deg);
end