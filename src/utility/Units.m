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
        minutes = 60;
        hours = 3600;
        days = 86400;

        radians = 1;
        degrees = pi / 180;

        kilograms = 1;
    end

    methods(Access = public)
        function [obj] = Units()
            error("Cannot construct instance of class Units.");
        end
    end

end