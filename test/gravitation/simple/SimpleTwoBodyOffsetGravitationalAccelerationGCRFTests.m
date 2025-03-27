classdef SimpleTwoBodyOffsetGravitationalAccelerationGCRFTests < matlab.unittest.TestCase
    methods (Test)
        function TestAgainstSymbollic(testCase)   

            r_sat  = sym("r_sat", [3, 1], ["real"]);
            r_moon = sym("r_moon", [3, 1], ["real"]);
            mu     = sym("mu", [1, 1], ["real"]);

            r_ddot = -mu * (r_sat - r_moon) / norm(r_sat - r_moon)^3;

            r_ddot_fun = matlabFunction(r_ddot, "Vars", {r_sat, r_moon, mu});
          
            mu_  = 4902.800066 .* Units.kilometers.^3 ./ Units.seconds.^2;
            r_sat_ = [1e5; 1e6; 2e4];
            r_moon_ = [2e2; 1e1; -3e2];

            r_ddot_ = SimpleTwoBodyOffsetGravitationalAccelerationGCRF((r_sat_ - r_moon_), mu_);
            expected_r_ddot_ = r_ddot_fun(r_sat_, r_moon_, mu_);
            
            testCase.verifyEqual( ...
                r_ddot_, ...
                expected_r_ddot_, ...
                'RelTol', 1e-9);  
        end
    end
end