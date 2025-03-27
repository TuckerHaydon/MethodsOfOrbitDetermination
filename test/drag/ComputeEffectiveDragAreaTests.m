classdef ComputeEffectiveDragAreaTests < matlab.unittest.TestCase
    methods (Test)
        function TestFullArea(testCase)

            % Traveling in X-direction
            satellite_velocity_orientation_gcrf = [1; 0; 0];

            % Facet oriented in X-direction
            facet_orientation_gcrf = [1; 0; 0];

            facet_area = 5;

            ignore_backwards = [false];

            [effective_area] = ComputeEffectiveDragArea( ...
                satellite_velocity_orientation_gcrf, ...
                facet_orientation_gcrf, ...
                facet_area, ...
                ignore_backwards);

            expected_effective_area = 5;

            testCase.verifyEqual( ...
                effective_area, ...
                expected_effective_area, ...
                'AbsTol', 1e-9); 
        end

        function TestNoArea(testCase)

            % Traveling in X-direction
            satellite_velocity_orientation_gcrf = [1; 0; 0];

            % Facet oriented in Y-direction
            facet_orientation_gcrf = [0; 1; 0];

            facet_area = 5;

            ignore_backwards = [false];

            [effective_area] = ComputeEffectiveDragArea( ...
                satellite_velocity_orientation_gcrf, ...
                facet_orientation_gcrf, ...
                facet_area, ...
                ignore_backwards);

            expected_effective_area = 0;

            testCase.verifyEqual( ...
                effective_area, ...
                expected_effective_area, ...
                'AbsTol', 1e-9); 
        end

        function TestPartialArea(testCase)

            % Traveling in X-direction
            satellite_velocity_orientation_gcrf = [1; 0; 0];

            % Facet oriented in XY-direction
            facet_orientation_gcrf = [1; 1; 0];
            facet_orientation_gcrf = facet_orientation_gcrf ./ norm(facet_orientation_gcrf);

            facet_area = 5;

            ignore_backwards = [false];

            [effective_area] = ComputeEffectiveDragArea( ...
                satellite_velocity_orientation_gcrf, ...
                facet_orientation_gcrf, ...
                facet_area, ...
                ignore_backwards);

            expected_effective_area = 5 * cosd(45);

            testCase.verifyEqual( ...
                effective_area, ...
                expected_effective_area, ...
                'AbsTol', 1e-9); 
        end

        function TestBackwardsArea(testCase)

            % Traveling in X-direction
            satellite_velocity_orientation_gcrf = [1; 0; 0];

            % Facet oriented in negative x-direction
            facet_orientation_gcrf = [-1; 0; 0];
            facet_orientation_gcrf = facet_orientation_gcrf ./ norm(facet_orientation_gcrf);

            facet_area = 5;

            ignore_backwards = [false];

            [effective_area] = ComputeEffectiveDragArea( ...
                satellite_velocity_orientation_gcrf, ...
                facet_orientation_gcrf, ...
                facet_area, ...
                ignore_backwards);

            expected_effective_area = 0;

            testCase.verifyEqual( ...
                effective_area, ...
                expected_effective_area, ...
                'AbsTol', 1e-9); 
        end

        function TestBackwardsAreaIgnored(testCase)

            % Traveling in X-direction
            satellite_velocity_orientation_gcrf = [1; 0; 0];

            % Facet oriented in negative x-direction
            facet_orientation_gcrf = [-1; 0; 0];
            facet_orientation_gcrf = facet_orientation_gcrf ./ norm(facet_orientation_gcrf);

            facet_area = 5;

            ignore_backwards = [true];

            [effective_area] = ComputeEffectiveDragArea( ...
                satellite_velocity_orientation_gcrf, ...
                facet_orientation_gcrf, ...
                facet_area, ...
                ignore_backwards);

            expected_effective_area = 5;

            testCase.verifyEqual( ...
                effective_area, ...
                expected_effective_area, ...
                'AbsTol', 1e-9); 
        end

        function TestMultipleFacets(testCase)

            % Traveling in X-direction
            satellite_velocity_orientation_gcrf = [1; 0; 0];

            % Facet oriented in negative x-direction
            facet_orientation_gcrf = transpose([ ...
                +1, 0, 0;
                -1, 0, 0;
                0, +1, 0;
                0, -1, 0;
                0, 0, +1;
                0, 0, -1;
                -1, 1, 0;
            ]);
            facet_orientation_gcrf = facet_orientation_gcrf ./ vecnorm(facet_orientation_gcrf, 2, 1);

            facet_area = [6; 6; 8; 8; 12; 12; 15];

            ignore_backwards = [false; false; false; false; false; false; true];

            [effective_area] = ComputeEffectiveDragArea( ...
                satellite_velocity_orientation_gcrf, ...
                facet_orientation_gcrf, ...
                facet_area, ...
                ignore_backwards);

            expected_effective_area = 6 + 15 * cosd(45);

            testCase.verifyEqual( ...
                effective_area, ...
                expected_effective_area, ...
                'AbsTol', 1e-9); 
        end

        function TestMultipleVector(testCase)

            % Traveling in X-direction
            satellite_velocity_orientation_gcrf = transpose([ ...
                +1, 0, 0; ...
                -1, 0, 0;
            ]);

            % Facet oriented in negative x-direction
            facet_orientation_gcrf = transpose([ ...
                +1, 0, 0;
                -1, 0, 0;
                0, +1, 0;
                0, -1, 0;
                0, 0, +1;
                0, 0, -1;
                -1, 1, 0;
            ]);
            facet_orientation_gcrf = facet_orientation_gcrf ./ vecnorm(facet_orientation_gcrf, 2, 1);
            facet_orientation_gcrf = repmat(facet_orientation_gcrf, [1, 1, 2]);

            facet_area = transpose([ ...
                6, 6, 8, 8, 12, 12, 15;
                3, 3, 4, 4, 6, 6, 7.5;
            ]);

            ignore_backwards = transpose([ ...
                false, false, false, false, false, false, true;
                false, false, false, false, false, false, false;
            ]);

            [effective_area] = ComputeEffectiveDragArea( ...
                satellite_velocity_orientation_gcrf, ...
                facet_orientation_gcrf, ...
                facet_area, ...
                ignore_backwards);

            expected_effective_area = [
                6 + 15 * cosd(45);
                3 + 7.5 * cosd(45);
            ];

            testCase.verifyEqual( ...
                effective_area, ...
                expected_effective_area, ...
                'AbsTol', 1e-9); 
        end
    end
end