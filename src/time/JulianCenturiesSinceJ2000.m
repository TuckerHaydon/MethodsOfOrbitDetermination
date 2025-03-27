function [julian_centuries] = JulianCenturiesSinceJ2000( ...
        julian_date)
    % Computes the approximate number of Julian centuries from the J2000.0 epoch. This function is often invoked with
    % UTC or UT1, so it may not actually represent the number of centuries since the epoch because the epoch is defined
    % differently in these difference time frames. Still, this conversion is often used.
    %
    % Requires:
    % - julian_date:
    %   - (:, :) double.
    %   - The Julian date to convert.
    % 
    % Returns:
    % - julian_centuries:
    %   - (:, :) double.
    %   - The corresponding Julian centuries.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, Fifth Edition
    arguments(Input)
        julian_date(:, :) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        julian_centuries(:, :) double {mustBeReal, mustBeFinite}
    end

    % Vallado, 3-43
    julian_centuries = (julian_date - 2451545.0) ./ 36525;
end