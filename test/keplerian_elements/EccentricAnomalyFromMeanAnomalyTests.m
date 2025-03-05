classdef EccentricAnomalyFromMeanAnomalyTests < matlab.unittest.TestCase
    methods (Test)
        function test1(testCase) 
            % Test 100 different angles and eccentricities. Only elliptical orbits.
            E = linspace(0, 2*pi - 1e-5, 100);
            e = linspace(1e-5, 1 - 1e-5, 100);

            [E_grid, e_grid] = ndgrid(E, e);

            % Compute the mean anomaly
            M_grid = MeanAnomalyFromEccentricAnomaly(E_grid, e_grid);

            % Compute the eccentric anomaly
            expected_E_grid = E_grid;
            E_grid = zeros(size(expected_E_grid));
            for idx = 1:numel(E_grid)
                E_grid(idx) = EccentricAnomalyFromMeanAnomaly(M_grid(idx), e_grid(idx));
            end

            testCase.verifyEqual( ...
                E_grid, ...
                expected_E_grid, ...
                'AbsTol', 1e-8);        
        end
    end
end