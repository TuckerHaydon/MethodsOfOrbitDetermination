classdef AttitudeTNRToGCRFTests < matlab.unittest.TestCase
    methods (Test)
        function TestBasicCase1(testCase)   
            % Velocity orthogonal to position

            satellite_position_gcrf = [1; 0; 0];
            satellite_velocity_gcrf = [0; 1; 0];

            [C_body2gcrf] = AttitudeTNRToGCRF( ...
                satellite_position_gcrf, ...
                satellite_velocity_gcrf);

            position_vectors_body = eye(3);

            position_vectors_gcrf = C_body2gcrf * position_vectors_body;

            expected_position_vectors_gcrf = transpose([
                0, 1, 0;
                0, 0, 1;
                1, 0, 0;
            ]);

            testCase.verifyEqual( ...
                position_vectors_gcrf, ...
                expected_position_vectors_gcrf, ...
                'AbsTol', 1e-9); 
        end

        function TestBasicCase2(testCase)   
            % Velocity tilted a bit into y

            satellite_position_gcrf = [1; 0; 0];
            satellite_velocity_gcrf = [-0.1; 1; 0];

            [C_body2gcrf] = AttitudeTNRToGCRF( ...
                satellite_position_gcrf, ...
                satellite_velocity_gcrf);

            position_vectors_body = eye(3);

            position_vectors_gcrf = C_body2gcrf * position_vectors_body;

            expected_position_vectors_gcrf = transpose([
                0, 1, 0;
                0, 0, 1;
                1, 0, 0;
            ]);

            testCase.verifyEqual( ...
                position_vectors_gcrf, ...
                expected_position_vectors_gcrf, ...
                'AbsTol', 1e-9); 
        end

        function TestBasicCase3(testCase)   
            % Velocity tilted a bit into z

            satellite_position_gcrf = [1; 0; 0];
            satellite_velocity_gcrf = [0; cos(0.1); sin(0.1)];

            [C_body2gcrf] = AttitudeTNRToGCRF( ...
                satellite_position_gcrf, ...
                satellite_velocity_gcrf);

            position_vectors_body = eye(3);

            position_vectors_gcrf = C_body2gcrf * position_vectors_body;

            expected_position_vectors_gcrf = transpose([
                0, cos(0.1), sin(0.1);
                0, -sin(0.1), cos(0.1);
                1, 0, 0;
            ]);

            testCase.verifyEqual( ...
                position_vectors_gcrf, ...
                expected_position_vectors_gcrf, ...
                'AbsTol', 1e-9); 
        end

        function TestMultipleVectors(testCase)   
            % Velocity tilted a bit into z

            satellite_position_gcrf = transpose([
                1, 0, 0;
                1, 0, 0;
            ]);

            satellite_velocity_gcrf = transpose([ ...
                0, cos(0.1), sin(0.1); ...
                -0.1, 1, 0; ...
            ]);

            [C_body2gcrf] = AttitudeTNRToGCRF( ...
                satellite_position_gcrf, ...
                satellite_velocity_gcrf);

            position_vectors_body = eye(3);

            position_vectors_gcrf = pagemtimes(C_body2gcrf, position_vectors_body);

            expected_position_vectors_gcrf(:, :, 1) = transpose([
                0, cos(0.1), sin(0.1);
                0, -sin(0.1), cos(0.1);
                1, 0, 0;
            ]);

            expected_position_vectors_gcrf(:, :, 2) = transpose([
                0, 1, 0;
                0, 0, 1;
                1, 0, 0;
            ]);

            testCase.verifyEqual( ...
                position_vectors_gcrf, ...
                expected_position_vectors_gcrf, ...
                'AbsTol', 1e-9); 
        end
    end
end