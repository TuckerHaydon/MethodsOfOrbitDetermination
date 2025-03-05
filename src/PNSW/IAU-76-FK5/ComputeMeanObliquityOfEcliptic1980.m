function [mean_obliquity_of_ecliptic_1980_deg] = ComputeMeanObliquityOfEcliptic1980( ...
        julian_centuries_since_epoch)
    % Compute the mean obliquity of the ecliptic using the 1980 IAU nutation model.
    % 
    % Requires:
    % - julian_centuries_since_epoch:
    %   - (:, :) double.
    %   - The number of Julian centuries since epoch.
    % 
    % Returns:
    % - mean_obliquity_of_ecliptic_1980_deg:
    %   - (:, :) double.
    %   - The mean obliquity of the ecliptic.
    %   - Unit degrees.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th edition
    arguments(Input)
        julian_centuries_since_epoch(:, :) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        mean_obliquity_of_ecliptic_1980_deg(:, :) double {mustBeReal, mustBeFinite}
    end

    % Rename
    TTT = julian_centuries_since_epoch;

    % Vallado, 3-84
    % Mean obliquity of the ecliptic 1980
    % I think the arcseconds approach maintains better accuracy.
    epsilon_bar_1980_arcseconds = 84381.448 - 46.8150 .* TTT - 0.00059 .* TTT.^2 + 0.001813 .* TTT.^3;
    
    % Convert arcseconds to degrees.
    epsilon_bar_1980_deg = epsilon_bar_1980_arcseconds * Conversions.ARCSECONDS_TO_DEGREES;

    % Wrap
    epsilon_bar_1980_deg = wrapTo180(epsilon_bar_1980_deg);

    % Output
    mean_obliquity_of_ecliptic_1980_deg = epsilon_bar_1980_deg;
end