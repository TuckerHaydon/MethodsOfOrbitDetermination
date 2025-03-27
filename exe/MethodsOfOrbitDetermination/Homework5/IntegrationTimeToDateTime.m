function [utc_date_time] = IntegrationTimeToDateTime(integration_time)
    % Converts integration time as a modified UTC Julian date in seconds to a UTC datetime object.
    % 
    % Requires:
    % - integration_time:
    %   - (:, :) double.
    %   - The integration time in seconds.
    % 
    % Returns:
    % - utc_date_time:
    %   - (:, :) datetime
    %   - The UTC time as a datetime object w/ no time zone.

    modified_julian_date = integration_time .* Conversions.SECONDS_TO_DAYS;

    utc_date_time = datetime(modified_julian_date, 'ConvertFrom', 'modifiedjuliandate');
end