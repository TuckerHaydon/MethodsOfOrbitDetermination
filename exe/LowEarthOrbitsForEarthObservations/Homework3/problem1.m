%% Preamble
clc; clear all; close all;

% References: 
% - https://en.wikipedia.org/wiki/Position_of_the_Sun
% - https://lweb.cfa.harvard.edu/~jzhao/times.html

% - https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/MPG%20Book/Release/Chapter2-TimeScales.pdf
% According to NASA, GMST is measured from mean equinox of date, not the epoch equinox.

%% Constants
degrees_per_day       = 360;
hours_per_day         = 24;
degrees_per_hour      = degrees_per_day / hours_per_day;
degrees_per_minute    = degrees_per_hour / 60;
degrees_per_second    = degrees_per_minute / 60;
degrees_per_arcminute = 1/60;
degrees_per_arcsecond = degrees_per_arcminute / 60;

%% Given
earth_gravitational_parameter = 3.98600415e14 * (Units.meters)^3 / (Units.seconds)^2;
earth_equatorial_radius       = 6378136.3 * Units.meters;
earth_angular_velocity        = 7.292115e-5 * Units.radians / Units.seconds;

epoch_year_local    = 2024;
epoch_month_local   = 4;
epoch_day_local     = 8;
epoch_hour_local    = 13;
epoch_minute_local  = 37;

%% Convert local time to UTC time.
% Austin is typically UTC-6, however with daylight savings time it was UTC-5.
epoch_date_time = datetime( ...
    epoch_year_local, ...
    epoch_month_local, ...
    epoch_day_local, ...
    epoch_hour_local, ...
    epoch_minute_local, ...
    0, ...
    'TimeZone', 'America/Chicago');

% Note that offset = Local Time - GMST so GMST = Local Time - Offset
[local_utc_offset, offset_from_normal_due_to_dst] = tzoffset(epoch_date_time);
delta_GMST_hours = hours(local_utc_offset);
delta_GMST_deg   = delta_GMST_hours * degrees_per_hour; 

epoch_year_utc    = epoch_year_local;
epoch_month_utc   = epoch_month_local;
epoch_day_utc     = epoch_day_local;
epoch_hour_utc    = epoch_hour_local - delta_GMST_hours; 
epoch_minute_utc  = epoch_minute_local;

% Check epoch hasn't wrapped around
assert(0 <= epoch_hour_utc && epoch_hour_utc <= 23);

%% Express UTC time in JD.
epoch_utc_jd = juliandate(epoch_date_time);

% %% Compute angles between J2K equinox and MOD equinox
% % Vallado, 3-43
% num_julian_centuries_since_j2k = (epoch_utc_jd - 2451545.0) / 36525;
% T_TT = num_julian_centuries_since_j2k;
% 
% % Vallado, 3-91
% zeta  = 2306.2181 * degrees_per_arcsecond * T_TT + 0.30188 * T_TT^2 + 0.017998 * T_TT^3;
% theta = 2004.3109 * degrees_per_arcsecond * T_TT - 0.42665 * T_TT^2 - 0.041833 * T_TT^3;
% z     = 2306.2181 * degrees_per_arcsecond * T_TT + 1.09468 * T_TT^2 + 0.018203 * T_TT^3;

%% Convert JD UTC to GMST
% Assuming UT1 and TT are close
% epoch_TT_JD = epoch_utc_jd;
% epoch_UT1_JD = epoch_utc_jd;
% [GMST_deg] = GMSTFromJD( ...
%         epoch_TT_JD, ...
%         epoch_UT1_JD); % deg

epoch_ut1_jd = epoch_utc_jd;
GMST_deg     = siderealTime(epoch_ut1_jd);

%% Compute ecliptic latitude / longitude and right ascension / declination
% n = number of days of TT from J2000.0.
n = epoch_utc_jd - 2451545.0;

% L = mean longitude of the sun
L = 280.460 + 0.9856474 * n; % Unit degrees
L = wrapTo360(L);

% g = mean anomaly
g = 357.528 + 0.9856003 * n; % Unit degrees
g = wrapTo360(g);

% lambda = ecliptic longitude
lambda = L + 1.915 * sind(g) + 0.020 * sind(2 * g); % Unit degrees

% beta = ecliptic_latitude
beta = 0; % Unit degrees

% epsilon = obliquity of the ecliptic
epsilon = 23.439 - 0.0000004 * n; % Unit degrees

% alpha = right ascension 
alpha = atand(cosd(epsilon) * tand(lambda)); % Unit degrees

% alpha is in same quadrant as lambda
if abs(lambda - alpha) > 180
    fprintf("Wrapping alpha because [lambda, alpha] = [%0.3f, %0.3f]\n", ...
        lambda, alpha)
    alpha = alpha + 180;
    alpha = wrapTo360(alpha);
end

assert(abs(alpha - lambda) < 90);

% Alternate calculation route
% f = 180 / pi;
% t = tand(epsilon/2)^2;
% alpha2 = lambda - f*t*sind(2*lambda) + (f/2)*t^2*sind(4*lambda);

% delta = declination
delta = asind(sind(epsilon) * sind(lambda)); % Unit degrees

% Rename
sun_ecliptic_longitude_deg = lambda;
sun_ecliptic_latitude_deg  = beta;
sun_right_ascension_deg    = alpha;
sun_declination_deg        = delta;

%% Compute geographic latitude / longitude
sun_geographic_latitude_deg  = sun_declination_deg;
sun_geographic_longitude_deg = sun_right_ascension_deg - GMST_deg;

fprintf("Geographic lat/lon of sun: [%0.3f, %0.3f] degrees\n", ...
    sun_geographic_latitude_deg, sun_geographic_longitude_deg);

%% Compute the lunar coordinates
% https://aa.usno.navy.mil/data/ssconf
%                              Austin, TX                              
%           Location:  W 97°43'48.0", N30°17'24.0",   150m           
%              (Longitude referred to Greenwich meridian)              
% 
% 
% 
%                        2024 Apr 08 18:37:00.0  (UT1) 
% 
% Object     R.A.     Dec.    Dist.   Z.D. Az.  Elong.   Diam.     Mag.
% 
%            h   m     °  '    km       °    °      °    '   "    Illum.
% Moon       1 11.6  + 7 36  353934    23  183   N  0   33 45.0      0%
%                          ** Sun in total eclipse **

sun_moon_separation_deg = 0;

moon_diameter_arcminute = 33;
moon_diameter_arcsecond = 45;
moon_diameter_deg = ...
    moon_diameter_arcminute * degrees_per_arcminute + ...
    moon_diameter_arcsecond * degrees_per_arcsecond;

moon_right_ascension_hour   = 1;
moon_right_ascension_minute = 11.6;
moon_right_ascension_deg = ...
    moon_right_ascension_hour * degrees_per_hour + ...
    moon_right_ascension_minute * degrees_per_minute;

moon_declination_deg       = +7;
moon_declination_arcminute = 36;
moon_declination_deg = ...
    moon_declination_deg + ...
    moon_declination_arcminute * degrees_per_arcminute;

moon_geographic_latitude_deg  = moon_declination_deg;
moon_geographic_longitude_deg = moon_right_ascension_deg - GMST_deg;

%% Reporting
fprintf("USNO indicates an angular separation of %d deg.\n", sun_moon_separation_deg);
fprintf("USNO indicates a moon angular diameter of %0.1f deg.\n", moon_diameter_deg);
fprintf("\n");

fprintf("The angular separation in RA and Dec of Sun - Moon is: \n");
fprintf("RA:    %06.3f - %06.3f = %0.3f deg\n", ...
    sun_right_ascension_deg, ...
    moon_right_ascension_deg, ...
    sun_right_ascension_deg - moon_right_ascension_deg)
fprintf("Dec:   %06.3f - %06.3f = %0.3f deg\n", ...
    sun_declination_deg, ...
    moon_declination_deg, ...
    sun_declination_deg - moon_declination_deg);
fprintf("\n");

fprintf("The geographic coordinates of the Sun and Moon are: \n");
fprintf("Sun [lat, lon]: [%07.3f, %07.3f] deg\n", ...
    sun_geographic_latitude_deg , ...
    sun_geographic_longitude_deg);
fprintf("Sun [lat, lon]: [%07.3f, %07.3f] deg\n", ...
    moon_geographic_latitude_deg , ...
    moon_geographic_longitude_deg);
fprintf("\n");


%% Helper function
function [GMST_deg] = GMSTFromJD( ...
        epoch_TT_JD, ...
        epoch_UT1_JD)
    
    % Rename
    JD_TT = epoch_TT_JD;
    JD_UT = epoch_UT1_JD;

    % Constants
    hours_per_day    = 24;
    days_per_century = 36525;
    degrees_per_hour = 3600 / 240;

    % Time of previous midnight
    JD0 = round(JD_UT) - 0.5; 

    % Hours elapsed since previous midnight
    H = (JD_UT - JD0) * hours_per_day;   
    
    % Number of days between epoch and J2000
    D_TT = JD_TT - 2451545.0; 

    % Number of days between previous midnight and J2000
    D_UT = JD0 - 2451545.0; 
    
    % Number of centuries since J2000
    T = D_TT / days_per_century; 
    
    % Compute the GMST in hours
    GMST_hours = mod(6.697375 + 0.065707485828 * D_UT + 1.0027379 * H + 0.0854103 * T + 0.0000258 * T^2, 24);

    % Convert hours to degrees
    GMST_deg = GMST_hours * degrees_per_hour;
end