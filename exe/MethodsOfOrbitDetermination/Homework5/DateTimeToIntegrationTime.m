function [integration_time] = DateTimeToIntegrationTime(utc_date_time)
    % Converts a UTC date time to an integration time. Integration time is a modified Julian date in unit seconds.
    %
    % Requires:
    % - utc_date_time:
    %   - (:, :) datetime
    %   - The UTC time as a datetime object w/ no time zone.
    % 
    % Returns:
    % - integration_time:
    %   - (:, :) double.
    %   - The integration time in seconds.
    arguments(Input)
        utc_date_time(:, :) datetime
    end

    arguments(Output)
        integration_time(:, :) double {mustBeReal, mustBeFinite}
    end
    
    assert(all(isempty(utc_date_time.TimeZone), "all"), "TimeZone must be empty.");

    % [~, modified_julian_date] = JulianDate( ...
    %     utc_date_time.Year, ...
    %     utc_date_time.Month, ...
    %     utc_date_time.Day, ...
    %     utc_date_time.Hour, ...
    %     utc_date_time.Minute, ...
    %     utc_date_time.Second);

    modified_julian_date = mjuliandate(utc_date_time);

    integration_time = modified_julian_date .* Conversions.DAYS_TO_SECONDS;
end