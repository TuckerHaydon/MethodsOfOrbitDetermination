classdef AttitudeITRF2PEFTests < matlab.unittest.TestCase
    methods (Test)
        function TestVector(testCase)   

            polar_motion = transpose([0.001, 0.002]);
            [C_itrf2pef] = AttitudeITRF2PEF(polar_motion);

            expected_C_itrf2pef = [
                1,       0,     -0.001;
                0,       1,      0.002;
                0.001,   -0.002,      1;
            ];

            % Accuracy lowered because of orthonormalization
            testCase.verifyEqual( ...
                C_itrf2pef, ...
                expected_C_itrf2pef, ...
                'AbsTol', 1e-4); 
        end

        function TestMultipleVectors(testCase)   

            polar_motion = transpose([ ...
                0.001,  0.002;
               -0.005,  0.001;
                0.002, -0.003;

            ]);
            [C_itrf2pef] = AttitudeITRF2PEF(polar_motion);

            expected_C_itrf2pef(:, :, 1) = [
                1,       0,     -0.001;
                0,       1,      0.002;
                0.001,   -0.002,      1;
            ];

            expected_C_itrf2pef(:, :, 2) = [
                1,       0,      0.005;
                0,       1,      0.001;
               -0.005,   -0.001,      1;
            ];

            expected_C_itrf2pef(:, :, 3) = [
                1,       0,     -0.002;
                0,       1,     -0.003;
                0.002,    0.003,      1;
            ];

            % Accuracy lowered because of orthonormalization
            testCase.verifyEqual( ...
                C_itrf2pef, ...
                expected_C_itrf2pef, ...
                'AbsTol', 1e-4); 
        end
    end
end