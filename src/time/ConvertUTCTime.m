function [converted_time] = ConvertUTCTime(...
        utc_time, ...
        ut1_utc_sec) 
    % Converts UTC time to the following other times:
    % (1) UTC (Universal coordinated time)
    % (2) UT1 (Universal time)
    % (3) TAI (International atomic time)
    % (4) TT (Terrestrial time)
    % (5) GPS (Global Positioning System) 
    % (6) TDB (barycentric dynamical time)
    % (7) TCB (Barycentric coordinate time)
    % (8) TCG (Geocentric coordinate time)
    %
    % All time values are provided as datetime objects.
    % 
    % Requires:
    % - utc_time:
    %   - (1, 1) datetime.
    %   - The current UTC time.
    % - ut1_utc_sec:
    %   - (1, 1) double.
    %   - The current UT1 offset from UTC.
    %   - Provided by Earth orientation parameters.
    % 
    % Returns:
    % - converted_time (struct):
    %   - utc:
    %     - (1, 1) datetime.
    %     - The UTC time.
    %   - ut1:
    %     - (1, 1) datetime.
    %     - The UT1 time.
    %   - tai:
    %     - (1, 1) datetime.
    %     - The TAI time.
    %   - tt:
    %     - (1, 1) datetime.
    %     - The TT time.
    %   - gps:
    %     - (1, 1) datetime.
    %     - The GPS time.
    %   - tdb:
    %     - (1, 1) datetime.
    %     - The TDB time.
    %   - tcb:
    %     - (1, 1) datetime.
    %     - The tcb time.
    %   - tcg:
    %     - (1, 1) datetime.
    %     - The tcg time.
    % 
    % References:
    % - Vallado, Fundamentals of Astrodynamics and Applications
    arguments(Input)
        utc_time(1, 1) datetime
        ut1_utc_sec(1, 1) double
    end

    arguments(Output)
        converted_time(1, 1) struct
    end

    % No time zone
    assert(isempty(utc_time.TimeZone));

    % Determine the number of leap seconds
    leapseconds_table     = leapseconds;
    leapseconds_table_idx = find(utc_time > leapseconds_table.Date, 1, 'last');
    leapseconds_at_epoch  = leapseconds_table.CumulativeAdjustment(leapseconds_table_idx);
    leapseconds_at_epoch  = ((2*(leapseconds_table.Type(leapseconds_table_idx) == "+")) - 1) * leapseconds_at_epoch;

    % Offset definitions
    tai_utc_sec = 10 + seconds(leapseconds_at_epoch); % 10 + leapseconds by definition
    tt_tai_sec  = 32.184; % fixed by definition
    gps_tai_sec = -19; % fixed by definition
   
    % Compute simple times
    ut1_time = utc_time + seconds(ut1_utc_sec);
    tai_time = utc_time + seconds(tai_utc_sec);
    tt_time  = tai_time + seconds(tt_tai_sec);
    gps_time = tai_time + seconds(gps_tai_sec);

    % Compute the more complex times
    tdb_time = BarycentricDynamicalTimeFromTerrestrialTime(tt_time);

    % Output
    converted_time.utc = utc_time;
    converted_time.ut1 = ut1_time;
    converted_time.tai = tai_time;
    converted_time.tt  = tt_time;
    converted_time.gps = gps_time;
    converted_time.tdb = tdb_time;
    % converted_time.tcb = tcb_time;
    % converted_time.tcg = tcg_time;
end