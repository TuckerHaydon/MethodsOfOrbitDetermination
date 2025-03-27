function [effective_area] = ComputeEffectiveDragArea( ...
        satellite_velocity_orientation_gcrf, ...
        facet_orientation_gcrf, ...
        facet_area, ...
        ignore_backwards)
    % Computes the effective area subject to drag on a satellite.
    %
    % Requires:
    % - satellite_velocity_orientation_gcrf:
    %   - (3, :) double.
    %   - The unit satellite velocity.
    % - facet_orientation_gcrf:
    %   - (3, :, :) double.
    %   - The unit facet orientation vectors.
    %   - Each page is associated with one velocity unit vector.
    % - facet_area:
    %   - (:, :) double.
    %   - The facet areas.
    %   - Unit m^2.
    %   - Each column is associated with one velocity unit vector.
    % - ignore_backwards:
    %   - (:, :) bool.
    %   - Whether to ignore whether a facet is 'backwards' in the area computation.
    %   - Facets with a negative angle relative to the velocity vector have their effective areas typically set to zero.
    %   - This flag indicates that they should not be set to zero and that both sides of the facet should be used in the
    %     area computation.
    %   - Each column is associated with one velocity unit vector.
    % 
    % Returns:
    % - effective_area:
    %   - (:, 1) double.
    %   - The sum total effective area.
    %   - Unit m^2.
    arguments(Input)
        satellite_velocity_orientation_gcrf(3, :) {mustBeReal, mustBeFinite}
        facet_orientation_gcrf(3, :, :) {mustBeReal, mustBeFinite}
        facet_area(:, :) {mustBeReal, mustBeFinite, mustBeNonnegative}
        ignore_backwards(:, :) logical
    end

    arguments(Output)
        effective_area(:, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    num_vectors = size(satellite_velocity_orientation_gcrf, 2);
    num_facets = size(facet_orientation_gcrf, 2);

    assert(num_vectors == size(facet_orientation_gcrf, 3));
    assert(num_vectors == size(facet_area, 2));
    assert(num_vectors == size(ignore_backwards, 2));

    assert(num_facets == size(facet_area, 1));
    assert(num_facets == size(ignore_backwards, 1));

    % Make sure the orientations are unit length.
    assert(all(abs(vecnorm(satellite_velocity_orientation_gcrf, 2, 1) - 1) < 1e-9, 'all')); 
    assert(all(abs(vecnorm(facet_orientation_gcrf, 2, 1) - 1) < 1e-9, 'all')); 

    % Compute the angle between each facet and the velocity vector
    % satellite_velocity_orientation_gcrf is [3, N]
    % facet_orientation_gcrf is [3, M, N]
    facet_angles = acos(pagemtimes(...
        permute(satellite_velocity_orientation_gcrf, [3, 1, 2]), ...
        facet_orientation_gcrf));
    facet_angles = permute(facet_angles, [2, 3, 1]);

    % The facets to be included in the area computation are those whose angle relative to the velocity vector is less
    % than 90 degrees. Or we're ignoring backwards angles.
    facet_mask = ignore_backwards | (facet_angles < pi/2);

    % The effective area is the area times the cosine of the incidence angle. At 0 degrees, the full area is exposed. At
    % 90 degrees, no area is exposed.
    effective_areas = facet_area .* cos(facet_angles);

    % Make sure all areas are positive. Sometimes they can be negative if backwards angles are included.
    effective_areas = abs(effective_areas);

    % Sum all the areas together.
    effective_area = 0 + transpose(sum(effective_areas .* facet_mask, 1));
end