classdef OrthonormalizeMatrixTests < matlab.unittest.TestCase
    methods (Test)
        function TestMatrix1(testCase)   

            C_in = eye(3) + diag([0.1, 0.1, 0.1]);
            C_out = OrthonormalizeMatrix(C_in);

            expected_C_out = eye(3);

            testCase.verifyEqual( ...
                C_out, ...
                expected_C_out, ...
                'AbsTol', 1e-9); 
        end

        function TestMatrix2(testCase)   

            C_in = eye(3) + [
                0.1, 0.01, 0.2;
                0.1, -0.01, 0.1;
                0.003, 0.2, 0.1;
            ];
            C_out = OrthonormalizeMatrix(C_in);

            testCase.verifyEqual( ...
                det(C_out), ...
                1, ...
                'AbsTol', 1e-9); 

            testCase.verifyEqual( ...
                C_out * transpose(C_out), ...
                eye(3), ...
                'AbsTol', 1e-9); 
        end

        function TestArray(testCase)   

            C_in(:, :, 1) = eye(3) + [0.01, 0.02, 0.01; 0, 0, 0; 0, 0, 0];
            C_in(:, :, 2) = eye(3) + [0, 0, 0; 0.01, 0.02, 0.01; 0, 0, 0];
            C_in(:, :, 3) = eye(3) + [0, 0, 0; 0.01, 0.02, 0.01; 0, 0, 0];

            C_out = OrthonormalizeMatrix(C_in);

            expected_C_out(:, :, 1) = OrthonormalizeMatrix(C_in(:, :, 1));
            expected_C_out(:, :, 2) = OrthonormalizeMatrix(C_in(:, :, 2));
            expected_C_out(:, :, 3) = OrthonormalizeMatrix(C_in(:, :, 3));

            testCase.verifyEqual( ...
                C_out, ...
                expected_C_out, ...
                'AbsTol', 1e-9); 
        end
    end
end