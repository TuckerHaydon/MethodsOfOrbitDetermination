classdef OrbitalPeriodFromSemimajorAxisTests < matlab.unittest.TestCase
    methods (Test)
        function test1(testCase)            
            mu = 1;
            a = 1;
            e = 0.5; % This doesnt do anything except check to make sure it's an orbit

            expected_T = 2*pi*sqrt(1 / 1);
            T = OrbitalPeriodFromSemimajorAxis(a, e, mu);

            testCase.verifyEqual( ...
                T, ...
                expected_T, ...
                'RelTol', 1e-6, ...
                'AbsTol', 1e-6);
        end

        function test2(testCase)            
            mu = 1;
            a = 2;
            e = 0.5; % This doesnt do anything except check to make sure it's an orbit

            expected_T = 2*pi*sqrt(8 / 1);
            T = OrbitalPeriodFromSemimajorAxis(a, e, mu);

            testCase.verifyEqual( ...
                T, ...
                expected_T, ...
                'RelTol', 1e-6, ...
                'AbsTol', 1e-6);
        end
    end
end