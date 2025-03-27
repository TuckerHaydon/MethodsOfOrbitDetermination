classdef GMSTFromUT1Tests < matlab.unittest.TestCase
    methods (Test)
        function TestVallado(testCase) 
            % Vallado example 3-5
            ut1_jd = 2448855.009722;

            num_julian_centuries_of_ut1_time_since_j2000 = ...
                (ut1_jd - 2451545.0) / 36525;

            GMST_deg = GMSTFromUT1(num_julian_centuries_of_ut1_time_since_j2000);
            expected_GMST_deg = 152.578787810;

            testCase.verifyEqual( ...
                GMST_deg, ...
                expected_GMST_deg, ...
                'AbsTol', 1e-4); 
        end

        function TestAgainstBuiltIn(testCase) 
            % Vallado example 3-5
            ut1_jd = 2448855.009722;

            num_julian_centuries_of_ut1_time_since_j2000 = ...
                (ut1_jd - 2451545.0) / 36525;

            GMST_deg = GMSTFromUT1(num_julian_centuries_of_ut1_time_since_j2000);
            
            utc_jd = ut1_jd;
            expected_GMST_deg = siderealTime(utc_jd);

            testCase.verifyEqual( ...
                GMST_deg, ...
                expected_GMST_deg, ...
                'AbsTol', 1e-4); 
        end
    end
end