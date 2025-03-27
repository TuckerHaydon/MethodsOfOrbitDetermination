classdef ComputeMoonPositionGCRFTests < matlab.unittest.TestCase
    methods (Test)
        function TestVallado(testCase)
            % Vallado, Example 5-3
            utc_time = datetime(1994, 04, 27, 23, 58, 59.816);
            ut1_utc = -0.088972;

            time_conversions = ConvertUTCTime(utc_time, ut1_utc);

            earth_radius = 6378136.3;
            [moon_position_gcrf] = ComputeMoonPositionGCRF( ...
                time_conversions.tdb, ...
                earth_radius);

            expected_moon_position_gcrf = [
                -134240.611;
                -311571.556;
                -126693.771;
            ] .* Units.kilometers;

            testCase.verifyEqual( ...
                moon_position_gcrf, ...
                expected_moon_position_gcrf, ...
                'RelTol', 1e-6, ...
                'AbsTol', 1e1); 
        end

    end
end