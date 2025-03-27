function [bool] = IsParabolicOrbit( ...
        eccentricity, ...
        options)
    % Determines if an orbit is parabolic from its eccentricity.
    %
    % Requires: 
    % - eccentricity:
    %   - (:, :) double.
    %   - The orbital eccentricity.
    % - options (struct):
    %   - tolerance:
    %     - (1, 1) double.
    %     - The eccentricity tolerance to within a parabola is determined.
    arguments(Input)
        eccentricity(:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        options.tolerance(1, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 1e-6
    end

    arguments(Output)
        bool(:, :) logical
    end

    bool = abs(eccentricity - 1.0) <= options.tolerance;

    % Output check
    assert(all(size(bool) == size(eccentricity)));
end