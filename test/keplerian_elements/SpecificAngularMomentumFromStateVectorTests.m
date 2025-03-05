classdef SpecificAngularMomentumFromStateVectorTests < matlab.unittest.TestCase
    methods (Test)
        function test1(testCase)
            % Test from Vallado orbital elements example
            position_eci = [6524.834; 6862.875; 6448.296] .* Units.kilometers;
            velocity_eci = [4.901327; 5.533756; -1.976341] .* Units.kilometers ./ Units.seconds;

            expected_specific_angular_momentum = [-49246.7; 44500.5; 2469.6] .* (Units.kilometers).^2 ./ Units.seconds;

            [specific_angular_momentum] = SpecificAngularMomentumFromStateVector( ...
                position_eci, ...
                velocity_eci);

            testCase.verifyEqual( ...
                specific_angular_momentum, ...
                expected_specific_angular_momentum, ...
                'RelTol', 1e-3);
        end
    end
end