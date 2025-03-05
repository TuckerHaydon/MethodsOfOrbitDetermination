classdef ECIJ2GravitationTests < matlab.unittest.TestCase
    methods (Test)
        function TestAgainstSymbollicJacobian(testCase)   
            %% Setup
            mu = sym("mu", ["real", "positive"]);
            J2 = sym("J2", ["real", "positive"]);
            Re = sym("Re", ["real", "positive"]);
            x = sym("x", [3, 1], "real");
            
            r = norm(x);
            sin_phi = x(3) / r;
            
            % Gravitational potential
            U = (mu / r) * (1 - J2 * (Re/r)^2 * ((3/2) * sin_phi^2 - 1/2));
            
            %% Compute the Jacobian
            jacobian_U = jacobian(U, x);
            jacobian_U_fun = matlabFunction(jacobian_U, 'Vars', transpose([mu; J2; Re; x]));
            
            %% Compare against hand-calculated jacobian with Monte Carlo
            mu_ = 398600.4*10^9;
            J2_ = 0.00108248;
            Re_ = 6378.145e3;
            
            num_iter = 1e5;

            % Calculate random position vectors
            position_eci = zeros(3, num_iter);
            for idx = 1:num_iter
                % https://math.libretexts.org/Courses/Mount_Royal_University/MATH_2200%3A_Calculus_for_Scientists_II/7%3A_Vector_Spaces/5.7%3A_Cylindrical_and_Spherical_Coordinates
                theta = 2 * pi * rand();
                phi   = pi * rand();
                r     = Re_;
            
                position_eci(:, idx) = [
                    r * sin(phi) * cos(theta);
                    r * sin(phi) * sin(theta);
                    r * cos(phi);
                ];
            end

            % Compute the gravitational acceleration vector with the function
            [gravitation_eci] = ECIJ2Gravitation( ...
                position_eci, ...
                mu_, ...
                J2_, ...
                Re_);

            % Compute the gravitational acceleration vector with the symbollic math
            analytical_gravitation_eci = transpose(jacobian_U_fun( ...
                mu_, ...
                J2_, ...
                Re_, ...
                transpose(position_eci(1, :)), ...
                transpose(position_eci(2, :)), ...
                transpose(position_eci(3, :))));

            differences = gravitation_eci - analytical_gravitation_eci;

            testCase.verifyTrue(all(abs(differences) < 1e-12, 'all')); 
        end
    end
end