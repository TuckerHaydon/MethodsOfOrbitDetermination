classdef SpecificMechanicalEnergyFromStateVectorTests < matlab.unittest.TestCase
    methods (Test)
        function test1(testCase)            
            r = [0; 0; 1];
            v = [1; 0; 0];
            mu = 1;

            expected_xi = 0.5 * 1^2 - 1 / 1;
            xi = SpecificMechanicalEnergyFromStateVector(r, v, mu);

            testCase.verifyEqual( ...
                xi, ...
                expected_xi, ...
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