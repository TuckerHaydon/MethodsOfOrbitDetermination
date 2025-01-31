function [bool] = IsEquatorialOrbit( ...
        inclination, ...
        options)
    % Determines if an orbit is circular from its eccentricity.
    %
    % Requires: 
    % - inclination:
    %   - (:, :) double.
    %   - The inclination of the orbit.
    %   - Unit radians
    % - options (struct):
    %   - tolerance:
    %     - (1, 1) double.
    %     - The eccentricity tolerance to within an equatorial orbit is determined.
    arguments(Input)
        inclination(:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        options.tolerance(1, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 1e-6
    end

    arguments(Output)
        bool(:, :) logical
    end

    bool = (inclination <= options.tolerance);

    % Output check
    assert(all(size(bool) == size(inclination)));
end