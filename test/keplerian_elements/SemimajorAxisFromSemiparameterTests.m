classdef SemimajorAxisFromSemiparameterTests < matlab.unittest.TestCase
    methods (Test)
        function test1(testCase)
            semiparameter = 11067.798 .* Units.kilometers;
            eccentricity  = 0.832853;

            expected_semimajor_axis = semiparameter ./ (1.0 - eccentricity.^2);


            [semimajor_axis] = SemimajorAxisFromSemiparameter( ...
                semiparameter, ...
                eccentricity);

            testCase.verifyEqual( ...
                semimajor_axis, ...
                expected_semimajor_axis, ...
                'RelTol', 1e-12, ...
                'AbsTol', 1e-12);
        end
    end
end