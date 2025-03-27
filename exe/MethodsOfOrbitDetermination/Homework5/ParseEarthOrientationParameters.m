function [...      
        ut1_utc_sec, ...
        polar_motion_deg, ...
        length_of_day_earth_orientation_parameter_sec, ...
        longitude_earth_orientation_parameter_deg, ...
        obliquity_earth_orientation_parameter_deg] = ...
    ParseEarthOrientationParameters(...
        eop_file, ...
        utc_time)
    % Parses the Earth orientation parameters (EOP) associated with a single time from a provided EOP file.
    % 
    % Requires:
    % - eop_file:
    %   - (1, 1) string.
    %   - Path the the EOP file.
    % - utc_time:
    %   - (1, 1) datetime.
    %   - The UTC date time of interest.
    %
    % Returns:
    % - ut1_utc_sec:
    %   - (1, 1) double.
    %   - The UT1 offset from UTC.
    %   - Unit seconds.
    % - polar_motion_deg:
    %   - (2, 1) double.
    %   - The polar motion parameter.
    %   - Unit degrees.
    % - length_of_day_earth_orientation_parameter_sec:
    %   - (1, 1) double
    %   - The length of day parameter.
    %   - Unit seconds.
    % - longitude_earth_orientation_parameter_deg:
    %   - (1, 1) double
    %   - The nutation in longitude parameter.
    %   - Unit degrees.
    % - obliquity_earth_orientation_parameter_deg:
    %   - (1, 1) double.
    %   - The nutation in obliquity parameter.
    %   - Unit degrees.
    % 
    % References:
    % - https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
    arguments(Input)
        eop_file(1, 1) string
        utc_time(1, 1) datetime
    end

    arguments(Output)
        ut1_utc_sec(1, 1) double {mustBeReal, mustBeFinite}
        polar_motion_deg(2, 1) double {mustBeReal, mustBeFinite}
        length_of_day_earth_orientation_parameter_sec(1, 1) double {mustBeReal, mustBeFinite}
        longitude_earth_orientation_parameter_deg(1, 1) double {mustBeReal, mustBeFinite}
        obliquity_earth_orientation_parameter_deg(1, 1) double {mustBeReal, mustBeFinite}
    end

    %% Determine the Earth orientation parameters of the provided date
    % Downloaded the parameters from https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
    % There are 185 characters in a full line.
    
    % Read all lines
    eop_lines = readlines(eop_file);
    
    % Strip trailing whitespace
    eop_lines = strip(eop_lines, 'right');
    
    % Remove any incomplete lines
    num_expected_characters = 185;
    eop_lines = eop_lines(strlength(eop_lines) == num_expected_characters);
    
    % Turn it into a character array
    eop_lines = char(eop_lines);
    
    %% Parse the EOP
    % I believe error is the standard deviation.
    % Col.#    Format  Quantity
    % -------  ------  -------------------------------------------------------------
    % 1-2      I2      year (to get true calendar year, add 1900 for MJD<=51543 or add 2000 for MJD>=51544)
    % 3-4      I2      month number
    % 5-6      I2      day of month
    % 7        X       [blank]
    % 8-15     F8.2    fractional Modified Julian Date (MJD UTC)
    % 16       X       [blank]
    % 17       A1      IERS (I) or Prediction (P) flag for Bull. A polar motion values
    % 18       X       [blank]
    % 19-27    F9.6    Bull. A PM-x (sec. of arc)
    % 28-36    F9.6    error in PM-x (sec. of arc)
    % 37       X       [blank]
    % 38-46    F9.6    Bull. A PM-y (sec. of arc)
    % 47-55    F9.6    error in PM-y (sec. of arc)
    % 56-57    2X      [blanks]
    % 58       A1      IERS (I) or Prediction (P) flag for Bull. A UT1-UTC values
    % 59-68    F10.7   Bull. A UT1-UTC (sec. of time)
    % 69-78    F10.7   error in UT1-UTC (sec. of time)
    % 79       X       [blank]
    % 80-86    F7.4    Bull. A LOD (msec. of time) -- NOT ALWAYS FILLED
    % 87-93    F7.4    error in LOD (msec. of time) -- NOT ALWAYS FILLED
    % 94-95    2X      [blanks]
    % 96       A1      IERS (I) or Prediction (P) flag for Bull. A nutation values
    % 97       X       [blank]
    % 98-106   F9.3    Bull. A dPSI (msec. of arc)
    % 107-115  F9.3    error in dPSI (msec. of arc)
    % 116      X       [blank]
    % 117-125  F9.3    Bull. A dEPSILON (msec. of arc)
    % 126-134  F9.3    error in dEPSILON (msec. of arc)
    % 135-144  F10.6   Bull. B PM-x (sec. of arc)
    % 145-154  F10.6   Bull. B PM-y (sec. of arc)
    % 155-165  F11.7   Bull. B UT1-UTC (sec. of time)
    % 166-175  F10.3   Bull. B dPSI (msec. of arc)
    % 176-185  F10.3   Bull. B dEPSILON (msec. of arc)
    eop.daily.year                                        = int32(str2num(eop_lines(:, 1:2)));
    eop.daily.month                                       = int32(str2num(eop_lines(:, 3:4)));
    eop.daily.day                                         = int32(str2num(eop_lines(:, 5:6)));
    % line 7 blank
    eop.daily.mjd_utc                                     = str2double(string(eop_lines(:, 8:15)));
    % line 16 blank
    eop.daily.polar_motion_type                           = eop_lines(:, 17);
    % line 18 blank
    eop.daily.polar_motion_x_arcseconds                   = str2double(string(eop_lines(:, 19:27)));
    eop.daily.polar_motion_x_arcseconds_error             = str2double(string(eop_lines(:, 28:36)));
    % line 37 blank
    eop.daily.polar_motion_y_arcseconds                   = str2double(string(eop_lines(:, 38:46)));
    eop.daily.polar_motion_y_arcseconds_error             = str2double(string(eop_lines(:, 47:55)));
    % lines 56, 57 blank
    eop.daily.ut1_utc_type                                = eop_lines(:, 58);
    eop.daily.ut1_utc_sec                                 = str2double(string(eop_lines(:, 59:68)));
    eop.daily.ut1_utc_sec_error                           = str2double(string(eop_lines(:, 69:78)));
    % line 79 blank
    eop.daily.lod_millisec                                = str2double(string(eop_lines(:, 80:86)));
    eop.daily.lod_millisec_error                          = str2double(string(eop_lines(:, 87:93)));
    % lines 94, 95 blank
    eop.daily.nutation_type                               = eop_lines(:, 96);
    % line 97 blank
    eop.daily.nutation_dPsi_milliarcsecond                = str2double(string(eop_lines(:, 98:106)));
    eop.daily.nutation_dPsi_milliarcsecond_error          = str2double(string(eop_lines(:, 107:115)));
    % line 116 blank
    eop.daily.nutation_dEpsilon_milliarcsecond            = str2double(string(eop_lines(:, 117:125)));
    eop.daily.nutation_dEpsilon_milliarcsecond_error      = str2double(string(eop_lines(:, 126:134)));
    eop.monthly.polar_motion_x_arcseconds                 = str2double(string(eop_lines(:, 135:144)));
    eop.monthly.polar_motion_y_arcseconds                 = str2double(string(eop_lines(:, 145:154)));
    eop.monthly.ut1_utc_sec                               = str2double(string(eop_lines(:, 155:165)));
    eop.monthly.nutation_dPsi_milliarcsecond              = str2double(string(eop_lines(:, 166:175)));
    eop.monthly.nutation_dEpsilon_milliarcsecond          = str2double(string(eop_lines(:, 176:185)));
    
    % Convert year into four digits
    mask_1900 = eop.daily.mjd_utc <= 51543;
    mask_2000 = ~mask_1900;
    eop.daily.year(mask_1900) = eop.daily.year(mask_1900) + 1900;
    eop.daily.year(mask_2000) = eop.daily.year(mask_2000) + 2000;
    
    %% Find the data of interest
    year_mask  = eop.daily.year  == year(utc_time);
    month_mask = eop.daily.month == month(utc_time);
    day_mask   = eop.daily.day   == day(utc_time);
    
    data_mask  = year_mask & month_mask & day_mask;
    
    assert(sum(data_mask) == 1, "More or less than 1 record found!");
    
    data_idx = find(data_mask, 1, 'first');
    
    polar_motion_arcseconds = [
        eop.daily.polar_motion_x_arcseconds(data_idx); 
        eop.daily.polar_motion_y_arcseconds(data_idx); 
    ];
    
    ut1_utc_sec           = eop.daily.ut1_utc_sec(data_idx);
    lod_ms                = eop.daily.lod_millisec(data_idx);
    dPsi_milli_arcseconds = eop.daily.nutation_dPsi_milliarcsecond(data_idx);
    dEps_milli_arcseconds = eop.daily.nutation_dEpsilon_milliarcsecond(data_idx);
    
    ut1_utc_sec                                     = ut1_utc_sec;
    polar_motion_deg                                = polar_motion_arcseconds * Conversions.ARCSECONDS_TO_DEGREES;
    length_of_day_earth_orientation_parameter_sec   = lod_ms / 1000;
    longitude_earth_orientation_parameter_deg       = dPsi_milli_arcseconds / 1000 * Conversions.ARCSECONDS_TO_DEGREES;
    obliquity_earth_orientation_parameter_deg       = dEps_milli_arcseconds / 1000 * Conversions.ARCSECONDS_TO_DEGREES;
end