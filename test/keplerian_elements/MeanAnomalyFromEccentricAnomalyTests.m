classdef MeanAnomalyFromEccentricAnomalyTests < matlab.unittest.TestCase
    methods (Test)
        function testScalar(testCase)            
            % Vallado, Example 2-1 
            E = 220.512074767522 .* Units.degrees;
            e = 0.4;

            expected_M = 235.4 .* Units.degrees;

            M = MeanAnomalyFromEccentricAnomaly(E, e);

            testCase.verifyEqual( ...
                M, ...
                expected_M, ...
                'RelTol', 1e-6, ...
                'AbsTol', 1e-6);
        end
    end
end