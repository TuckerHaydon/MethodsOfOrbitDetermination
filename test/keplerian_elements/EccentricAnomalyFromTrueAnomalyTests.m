classdef EccentricAnomalyFromTrueAnomalyTests < matlab.unittest.TestCase
    methods (Test)
        function test1(testCase)            
            eccentricity = 0.5;
            true_anomaly = deg2rad(70);

            expected_eccentric_anomaly = atan2(...
                (sin(true_anomaly) * sqrt(1 - eccentricity^2)) / (1 + eccentricity * cos(true_anomaly)), ...
                (eccentricity + cos(true_anomaly)) / (1 + eccentricity * cos(true_anomaly)));
    
            eccentric_anomaly = EccentricAnomalyFromTrueAnomaly( ...
                true_anomaly, ...
                eccentricity);

            testCase.verifyEqual( ...
                eccentric_anomaly, ...
                expected_eccentric_anomaly, ...
                'AbsTol', 1e-8); 
        end
    end
end