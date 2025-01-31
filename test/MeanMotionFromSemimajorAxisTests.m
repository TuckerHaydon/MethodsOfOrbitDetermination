classdef MeanMotionFromSemimajorAxisTests < matlab.unittest.TestCase
    methods (Test)
        function test1(testCase)            
            mu = 1;
            a = 1;
            e = 0.5; % This doesnt do anything except check to make sure it's an orbit

            expected_n = 1;

            n = MeanMotionFromSemimajorAxis(a, e, mu);

            testCase.verifyEqual( ...
                n, ...
                expected_n, ...
                'RelTol', 1e-6, ...
                'AbsTol', 1e-6);
        end

        function test2(testCase)            
            mu = 1;
            a = 2;
            e = 0.5; % This doesnt do anything except check to make sure it's an orbit

            expected_n = sqrt(1/8);

            n = MeanMotionFromSemimajorAxis(a, e, mu);

            testCase.verifyEqual( ...
                n, ...
                expected_n, ...
                'RelTol', 1e-6, ...
                'AbsTol', 1e-6);
        end
    end
end