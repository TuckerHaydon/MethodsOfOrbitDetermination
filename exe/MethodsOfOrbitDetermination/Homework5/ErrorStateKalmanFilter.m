classdef ErrorStateKalmanFilter < handle
    properties(GetAccess = public, SetAccess = public)
        time
        nominal_state
        error_covariance
        dynamics_functions
        error_state_dynamics_functions
    end

    methods(Access = public)
        function [obj] = ErrorStateKalmanFilter( ...
                initial_time, ...
                initial_state, ...
                initial_error_covariance)

            obj.time             = initial_time;
            obj.nominal_state    = initial_state;
            obj.error_covariance = initial_error_covariance;

        end

        function [] = AttachSuperimposedDynamicsFunction( ...
                dynamics_function)
            % Attaches a dynamics function to the Kalman filter
        end
    end
end