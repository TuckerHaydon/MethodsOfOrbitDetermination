classdef CubatureTransformTests < matlab.unittest.TestCase
    methods (Test)
        function TestLinearSystem1(testCase) 
            % Linear transform should be exact.

            H_true = [
                4, -1;
            ];

            x = [
                1; 
                2;
            ];

            Pxx = [
                4, 1;
                1, 6;
            ];

            h = @(x) H_true * x;
            
            R = [
                    1
                ];

            expected_y = H_true * x;
            expected_Pyy = H_true * Pxx * transpose(H_true) + R;
            expected_Pxy = Pxx * transpose(H_true);
            expected_H = H_true;

            [y, Pyy, Pxy, H] = CubatureTransform(x, Pxx, h, R);
            
            testCase.verifyEqual( ...
                y, ...
                expected_y, ...
                'AbsTol', 1e-9); 

            testCase.verifyEqual( ...
                Pyy, ...
                expected_Pyy, ...
                'AbsTol', 1e-9); 

            testCase.verifyEqual( ...
                Pxy, ...
                expected_Pxy, ...
                'AbsTol', 1e-9); 

            testCase.verifyEqual( ...
                H, ...
                expected_H, ...
                'AbsTol', 1e-9); 
        end

        function TestLinearSystem2(testCase) 
            % Linear transform should be exact.

            H_true = [
                4, -1;
                2, 3;
                6, 2;
            ];

            x = [
                0; 
                5;
            ];

            Pxx = [
                4, 2;
                2, 6;
            ];

            h = @(x) H_true * x;
            
            R = [
                    1, 0.5, 0.1;
                    0.5, 1, 0.5;
                    0.1, 0.5, 2;
                ];

            expected_y = H_true * x;
            expected_Pyy = H_true * Pxx * transpose(H_true) + R;
            expected_Pxy = Pxx * transpose(H_true);
            expected_H = H_true;

            [y, Pyy, Pxy, H] = CubatureTransform(x, Pxx, h, R);
            
            testCase.verifyEqual( ...
                y, ...
                expected_y, ...
                'AbsTol', 1e-9); 

            testCase.verifyEqual( ...
                Pyy, ...
                expected_Pyy, ...
                'AbsTol', 1e-9); 

            testCase.verifyEqual( ...
                Pxy, ...
                expected_Pxy, ...
                'AbsTol', 1e-9); 

            testCase.verifyEqual( ...
                H, ...
                expected_H, ...
                'AbsTol', 1e-9); 
        end

        function TestMildNonlinearSystem(testCase) 
            % Mild nonlinear range case

            % Point at origin
            x = [0; 0];

            % Point far away
            u = [2000; 1000];

            h = @(x) vecnorm(u - x, 2, 1);

            Pxx = [
                4, 2;
                2, 6;
            ];

            R = [
                    1
                ];

            % Analytical H
            expected_H   = -transpose(u - x) ./ norm(u - x);

            expected_y   = h(x);
            expected_Pyy = expected_H * Pxx * transpose(expected_H) + R;
            expected_Pxy = Pxx * transpose(expected_H);

            [y, Pyy, Pxy, H] = CubatureTransform(x, Pxx, h, R);
            
            testCase.verifyEqual( ...
                y, ...
                expected_y, ...
                'AbsTol', 1e-3); 

            testCase.verifyEqual( ...
                Pyy, ...
                expected_Pyy, ...
                'AbsTol', 1e-3); 

            testCase.verifyEqual( ...
                Pxy, ...
                expected_Pxy, ...
                'AbsTol', 1e-3); 

            testCase.verifyEqual( ...
                H, ...
                expected_H, ...
                'AbsTol', 1e-3); 
        end
    end
end