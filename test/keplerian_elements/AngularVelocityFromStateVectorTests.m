classdef AngularVelocityFromStateVectorTests < matlab.unittest.TestCase
    methods (Test)
        function testScalar(testCase)   

            r = [3; 0; 0];
            v = [0; 2; 0];

            expected_angular_velocity = [0; 0; 6/9];

            angular_velocity = AngularVelocityFromStateVector(r, v);

            testCase.verifyEqual( ...
                angular_velocity, ...
                expected_angular_velocity, ...
                'AbsTol', 1e-8); 
        end

        function testVector(testCase)   

            r = transpose([
                3, 0, 0;
                0, 4, 0;
            ]);

            v = transpose([
                0, 2, 0;
                0, 0, 2;
            ]);

            expected_angular_velocity = transpose([
                0, 0, 6/9;
                8/16, 0, 0;
            ]);

            angular_velocity = AngularVelocityFromStateVector(r, v);

            testCase.verifyEqual( ...
                angular_velocity, ...
                expected_angular_velocity, ...
                'AbsTol', 1e-8); 
        end
    end
end