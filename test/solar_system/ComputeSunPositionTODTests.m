classdef ComputeSunPositionTODTests < matlab.unittest.TestCase
    methods (Test)
        function TestVallado(testCase)
            % Vallado, Example 5-1
            utc_time = datetime(2006, 04, 01, 23, 58, 54.816);
            ut1_utc = 0.265363;

            time_conversions = ConvertUTCTime(utc_time, ut1_utc);

            % The example approximates tdb with tt
            [sun_poisition_tod] = ComputeSunPositionTOD( ...
                time_conversions.ut1, ...
                time_conversions.tdb);

            expected_sun_position_tod = [
                146186588;
                28787227;
                12480305;
            ] .* Units.kilometers;

            testCase.verifyEqual( ...
                sun_poisition_tod, ...
                expected_sun_position_tod, ...
                'RelTol', 1e-4); 
        end

    end
end