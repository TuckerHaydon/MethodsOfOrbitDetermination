classdef Constants
    properties (Constant)
        SPEED_OF_LIGHT = 299792458;
    end

    methods(Access = public)
        function [obj] = Constants()
            error("Cannot construct instance of class Constants.");
        end
    end

end