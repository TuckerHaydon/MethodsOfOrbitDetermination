function [bool] = IsEllipticalOrbit( ...
        eccentricity, ...
        options)
    % Determines if an orbit is elliptical from its eccentricity.
    %
    % Requires: 
    % - eccentricity:
    %   - (:, :) double.
    %   - The orbital eccentricity.
    % - options (struct):
    %   - circle_tolerance:
    %     - (1, 1) double.
    %     - The eccentricity tolerance to within a parabola is determined.
    %   - parabola_tolerance:
    %     - (1, 1) double.
    %     - The eccentricity tolerance to within a parabola is determined.
    arguments(Input)
        eccentricity(:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        options.circle_tolerance(1, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 1e-6
        options.parabola_tolerance(1, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 1e-6
    end

    arguments(Output)
        bool(:, :) logical
    end

    bool = ...
        (eccentricity >= options.circle_tolerance) & ...
        (eccentricity <= 1 - options.parabola_tolerance);

    % Output check
    assert(all(size(bool) == size(eccentricity)));
end