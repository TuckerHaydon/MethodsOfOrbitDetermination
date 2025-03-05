function [GAST_deg] = GASTFromGMST( ...
        GMST_deg, ...
        equation_of_equinoxes_deg)
    % Computes the Greenwich apparent sidereal time (GAST) from the Greenwich mean sidereal time (GMST).
    %
    % Requires:
    % - GMST_deg:
    %   - (:, :) double.
    %   - The GMST.
    %   - Unit degrees.
    % - equation_of_equinoxes_deg:
    %   - (:, :) double.
    %   - The equation of the equinoxes.
    %   - Unit degrees.
    % 
    % Returns:
    % - GAST_deg:
    %   - (:, :) double.
    %   - The GAST.
    %   - Unit degrees.
    % 
    % References:
    % - Vallado, Fundamentals of astrodynamics and applications, 5th edition
    arguments(Input)
        GMST_deg(:, :) double {mustBeReal, mustBeFinite}
        equation_of_equinoxes_deg(:, :) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        GAST_deg(:, :) double {mustBeReal, mustBeFinite}
    end

    assert(all(size(GMST_deg) == size(equation_of_equinoxes_deg)), ...
        "Input must be same size!");
    
    % Vallado, 3-82
    GAST_deg = GMST_deg + equation_of_equinoxes_deg;
    
    % Wrap to 360
    GAST_deg = wrapTo360(GAST_deg);
end