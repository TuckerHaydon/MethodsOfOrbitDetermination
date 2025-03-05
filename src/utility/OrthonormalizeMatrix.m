function [C_out] = OrthonormalizeMatrix(C_in, num_iterations)
    % Renormalizes an orthonormal matrix.
    % 
    % Requires:
    % - C_in:
    %   - (3, 3, :) double.
    %   - The matrix to be renormalized.
    % - num_iterations:
    %   - (1, 1) double.
    %   - The number of iterations to use while orthonormalizing the matrix.
    %   - Usually 4 is more than sufficient for 16 digits of precision.
    % 
    % Returns:
    % - C_out:
    %   - (3, 3, :) double.
    %   - The renormalized matrix.
    % 
    % References:
    % - An Iterative Algorithm for Computing the Best Estimate of an Orthogonal Matrix
    % - https://www.jstor.org/stable/2949484?seq=1
    arguments(Input)
        C_in(3, 3, :) double {mustBeReal, mustBeFinite}
        num_iterations(1, 1) double {mustBeReal, mustBeFinite, mustBePositive} = 4
    end

    arguments(Output)
        C_out(3, 3, :) double {mustBeReal, mustBeFinite}
    end

    % Rename
    R = C_in;

    % Bjorck, section 3.
    for idx = 1:num_iterations
        R = 1.5 * R - 0.5 * pagemtimes(R, pagemtimes(R, 'transpose', R, 'none'));
    end

    % Output
    C_out = R;
end