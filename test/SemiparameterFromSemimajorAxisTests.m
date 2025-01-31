classdef SemiparameterFromSemimajorAxisTests < matlab.unittest.TestCase
    methods (Test)
        function test1(testCase)
            semimajor_axis = 36127.338 .* Units.kilometers;
            eccentricity   = 0.832853;

            expected_semiparameter = semimajor_axis .* (1.0 - eccentricity.^2);

            [semiparameter] = SemiparameterFromSemimajorAxis( ...
                semimajor_axis, ...
                eccentricity);

            testCase.verifyEqual( ...
                semiparameter, ...
                expected_semiparameter, ...
                'RelTol', 1e-12, ...
                'AbsTol', 1e-12);
        end
    end
end