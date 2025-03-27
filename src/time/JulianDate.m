function [julian_date, modified_julian_date] = JulianDate( ...
        year, ...
        month, ...
        day, ...
        hour, ...
        minute, ...
        second) 
    % Computes the Julian date and modified Julian date.
    %
    % Note, this function only works between 1900 and 2100.
    % 
    % Requires:
    % - year:
    %   - (:, :) double.
    %   - The year.
    % - month:
    %   - (:, :) double.
    %   - The month.
    % - day:
    %   - (:, :) double.
    %   - The day.
    % - hour:
    %   - (:, :) double.
    %   - The hour.
    % - minute:
    %   - (:, :) double.
    %   - The minute.
    % - second:
    %   - (:, :) double.
    %   - The second.
    % 
    % Returns:
    % - julian_date:
    %   - (:, :) double.
    %   - The julian date.
    % - modified_julian_date:
    %   - (:, :) double.
    %   - The modified julian date.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications, 5th edition
    arguments(Input)
        year(:, :) double {mustBeReal, mustBeFinite}
        month(:, :) double {mustBeReal, mustBeFinite}
        day(:, :) double {mustBeReal, mustBeFinite}
        hour(:, :) double {mustBeReal, mustBeFinite}
        minute(:, :) double {mustBeReal, mustBeFinite}
        second(:, :) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        julian_date(:, :) double {mustBeReal, mustBeFinite}
        modified_julian_date(:, :) double {mustBeReal, mustBeFinite}
    end

    assert(all(size(year) == size(month)));
    assert(all(size(year) == size(day)));
    assert(all(size(year) == size(hour)));
    assert(all(size(year) == size(minute)));
    assert(all(size(year) == size(second)));

    % The leveraged algorithm only works between the years 1900 and 2100.
    assert(...
        all(year > 1900 & year < 2100), ...
        "Function domain only defined between 1900 and 2100.");

    % Vallado, Algorithm 14
    julian_date = 367 .* year ...
        - floor((7 .* (year + floor((month + 9) ./ 12))) ./ 4) ...
        + floor((275 .* month) ./ 9) ...
        + day ...
        + 1721013.5 + ...
        (3600 .* hour + 60 .* minute + second) ./ 86400;

    modified_julian_date = julian_date - 2400000.5;
end