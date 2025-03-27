function [C_tnr2gcrf] = AttitudeTNRToGCRF( ...
        satellite_position_gcrf, ...
        satellite_velocity_gcrf)
    % Computes the frame transform matrix that converts a vector expressed in the Transverse-Normal-Radial (TNR) frame
    % to a vector expressed in the GCRF frame. The TNR is defined as:
    %   - Positive Z pointed radially away from Earth
    %   - Positive X aligned with primary flight direction and orthogonal to Z
    %   - Positive Y orthogonal to X and Z and forming an orthonormal triad.
    % 
    % Requires:
    % - satellite_position_gcrf:
    %   - (3, 1) double.
    %   - The satellite position in the GCRF frame.
    %   - Unit meters.
    % - satellite_velocity_gcrf:
    %   - (3, 1) double.
    %   - The satellite velocity in the GCRF frame.
    %   - Unit meters / second.
    % 
    % Returns:
    % - C_tnr2gcrf:
    %   - (3, 3) double.
    %   - A frame transfomation matrix that transforms a vector expressed in the TNR frame to one expressed in the GCRF
    %     frame.
    arguments(Input)
        satellite_position_gcrf(3, :) double {mustBeReal, mustBeFinite}
        satellite_velocity_gcrf(3, :) double {mustBeReal, mustBeFinite}
    end

    arguments(Output)
        C_tnr2gcrf(3, 3, :) double {mustBeReal, mustBeFinite}
    end

    % Compute unit vectors
    r_hat = satellite_position_gcrf ./ vecnorm(satellite_position_gcrf, 2, 1);
    v_hat = satellite_velocity_gcrf ./ vecnorm(satellite_velocity_gcrf, 2, 1);

    % Compute the tnr unit vectors in the GCRF frame
    plus_z_tnr_gcrf = r_hat;
    plus_y_tnr_gcrf = cross(r_hat, v_hat, 1);
    plus_x_tnr_gcrf = cross(plus_y_tnr_gcrf, plus_z_tnr_gcrf, 1);

    % Normalize the vectors as they may no longer be unit length.
    plus_z_tnr_gcrf = plus_z_tnr_gcrf ./ vecnorm(plus_z_tnr_gcrf, 2, 1);
    plus_y_tnr_gcrf = plus_y_tnr_gcrf ./ vecnorm(plus_y_tnr_gcrf, 2, 1);
    plus_x_tnr_gcrf = plus_x_tnr_gcrf ./ vecnorm(plus_x_tnr_gcrf, 2, 1);
    
    % Compose the frame transform matrix
    C_tnr2gcrf = reshape([
        plus_x_tnr_gcrf;
        plus_y_tnr_gcrf;
        plus_z_tnr_gcrf;
    ], 3, 3, []);

end