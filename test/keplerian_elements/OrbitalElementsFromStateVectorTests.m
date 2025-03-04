classdef OrbitalElementsFromStateVectorTests < matlab.unittest.TestCase
    methods (Test)
        function test1(testCase)
            % Example 19.1 from https://www.spacerl.com/wp-content/uploads/2019/01/Orbital-Elements.pdf
        
            position_eci = [-6045; -3490; 2500] .* Units.kilometers;
            velocity_eci = [-3.457; 6.618; 2.533] .* Units.kilometers ./ Units.seconds;
            gravitational_parameter = 398600 .* (Units.kilometers).^3 ./ (Units.seconds).^2;
        
            expected_semimajor_axis                 = 8788 * Units.kilometers;
            expected_eccentricity                   = 0.1712;
            expected_inclination                    = 153.2 * Units.degrees;
            expected_longitude_of_ascending_node    = 255.3 * Units.degrees;
            expected_argument_of_periapsis          = 20.07 * Units.degrees;
            expected_true_anomaly                   = 28.45 * Units.degrees;

            [...
                ~, ...
                semimajor_axis, ...
                eccentricity, ...
                inclination, ...
                longitude_of_ascending_node, ...
                argument_of_periapsis, ...
                true_anomaly, ...
                ~, ...
                ~, ...
                ~] = ...
            OrbitalElementsFromStateVector( ...
                position_eci, ...
                velocity_eci, ...
                gravitational_parameter);

            testCase.verifyEqual( ...
                semimajor_axis, ...
                expected_semimajor_axis, ...
                'RelTol', 1e-2);
            testCase.verifyEqual( ...
                eccentricity, ...
                expected_eccentricity, ...
                'RelTol', 1e-2);
            testCase.verifyEqual( ...
                inclination, ...
                expected_inclination, ...
                'RelTol', 1e-2);
            testCase.verifyEqual( ...
                longitude_of_ascending_node, ...
                expected_longitude_of_ascending_node, ...
                'RelTol', 1e-2);
            testCase.verifyEqual( ...
                argument_of_periapsis, ...
                expected_argument_of_periapsis, ...
                'RelTol', 1e-2);
            testCase.verifyEqual( ...
                true_anomaly, ...
                expected_true_anomaly, ...
                'RelTol', 1e-2);
        end

        function test2(testCase)
            % From Vallado's textbook
        
            position_eci = [6524.834; 6862.875; 6448.297] .* Units.kilometers;
            velocity_eci = [4.901327; 5.533756; -1.976341] .* Units.kilometers ./ Units.seconds;
            gravitational_parameter = 398600 .* (Units.kilometers).^3 ./ (Units.seconds).^2;
        
            expected_semiparameter                  = 11067.798 * Units.kilometers;
            expected_semimajor_axis                 = 36127.338 * Units.kilometers;
            expected_eccentricity                   = 0.832853;
            expected_inclination                    = 87.870 * Units.degrees;
            expected_longitude_of_ascending_node    = 227.898 * Units.degrees;
            expected_argument_of_periapsis          = 53.38 * Units.degrees;
            expected_true_anomaly                   = 92.335 * Units.degrees;
            expected_argument_of_latitude           = 145.72009 * Units.degrees;
            expected_true_longitude                 = 55.282587 * Units.degrees;
            expected_true_longitude_of_periapsis    = 247.806 * Units.degrees; 

            [...
                semiparameter, ...
                semimajor_axis, ...
                eccentricity, ...
                inclination, ...
                longitude_of_ascending_node, ...
                argument_of_periapsis, ...
                true_anomaly, ...
                argument_of_latitude, ...
                true_longitude, ...
                true_longitude_of_periapsis] = ...
            OrbitalElementsFromStateVector( ...
                position_eci, ...
                velocity_eci, ...
                gravitational_parameter);

            testCase.verifyEqual( ...
                semiparameter, ...
                expected_semiparameter, ...
                'RelTol', 1e-4);
            testCase.verifyEqual( ...
                semimajor_axis, ...
                expected_semimajor_axis, ...
                'RelTol', 1e-4);
            testCase.verifyEqual( ...
                eccentricity, ...
                expected_eccentricity, ...
                'RelTol', 1e-4);
            testCase.verifyEqual( ...
                inclination, ...
                expected_inclination, ...
                'RelTol', 1e-4);
            testCase.verifyEqual( ...
                longitude_of_ascending_node, ...
                expected_longitude_of_ascending_node, ...
                'RelTol', 1e-4);
            testCase.verifyEqual( ...
                argument_of_periapsis, ...
                expected_argument_of_periapsis, ...
                'RelTol', 1e-4);
            testCase.verifyEqual( ...
                true_anomaly, ...
                expected_true_anomaly, ...
                'RelTol', 1e-4);
            testCase.verifyEqual( ...
                argument_of_latitude, ...
                expected_argument_of_latitude, ...
                'RelTol', 1e-4);
            testCase.verifyEqual( ...
                true_longitude, ...
                expected_true_longitude, ...
                'RelTol', 1e-4);
            testCase.verifyEqual( ...
                true_longitude_of_periapsis, ...
                expected_true_longitude_of_periapsis, ...
                'RelTol', 1e-4);
        end
    end
end