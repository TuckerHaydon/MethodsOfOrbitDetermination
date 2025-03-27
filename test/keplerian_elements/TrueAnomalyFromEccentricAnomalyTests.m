classdef TrueAnomalyFromEccentricAnomalyTests < matlab.unittest.TestCase
    methods (Test)
        function test1(testCase)            
            E = [pi];
            e = [0.5];

            sin_nu = sin(E) .* sqrt(1 - e.^2) / (1 - e .* cos(E));
            cos_nu = (cos(E) - e) / (1 - e .* cos(E));
            expected_nu = atan2(sin_nu, cos_nu);

            nu = TrueAnomalyFromEccentricAnomaly(E, e);

            testCase.verifyEqual( ...
                nu, ...
                expected_nu, ...
                'RelTol', 1e-6, ...
                'AbsTol', 1e-6);
        end

        function test2(testCase)            
            r = [0; 0; 4];
            v = [2; 0; 0];
            mu = 1;

            expected_xi = 0.5 * 2^2 - 1 / 4;
            xi = SpecificMechanicalEnergyFromStateVector(r, v, mu);

            testCase.verifyEqual( ...
                xi, ...
                expected_xi, ...
                'RelTol', 1e-6, ...
                'AbsTol', 1e-6);
        end
    end
end