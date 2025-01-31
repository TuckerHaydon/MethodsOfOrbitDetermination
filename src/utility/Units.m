classdef Units
    properties (Constant)
        hertz = 1;
        kilohertz = 1e3;
        megahertz = 1e6;
        gigahertz = 1e9;

        meters = 1;
        centimeters = 1e-2;
        millimeters = 1e-3;
        kilometers = 1e3;

        seconds = 1;

        radians = 1;
        degrees = pi / 180;
    end

    methods(Access = public)
        function [obj] = Units()
            error("Cannot construct instance of class Units.");
        end
    end

end