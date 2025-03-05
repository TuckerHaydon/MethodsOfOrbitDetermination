classdef Conversions
    properties (Constant)

        DEGREES_TO_HOURS = 24 / 360;
        HOURS_TO_DEGREES = 360 / 24;

        DEGREES_TO_MINUTES = (24 / 360) * (60 / 1);
        MINUTES_TO_DEGREES = (1 / 60) * (360 / 24);

        DEGREES_TO_SECONDS = (24 / 360) * (60 / 1) * (60 / 1);
        SECONDS_TO_DEGREES = (1 / 60) * (1 / 60) * (360 / 24);

        DEGREES_TO_ARCMINUTES = 60;
        ARCMINUTES_TO_DEGREES = 1 / 60;

        DEGREES_TO_ARCSECONDS = 60 * 60;
        ARCSECONDS_TO_DEGREES = 1 / (60 * 60);

        REVOLUTIONS_TO_DEGREES = 360;
        DEGREES_TO_REVOLUTIONS = 1 / 360;

        HOURS_TO_MINUTES = 60;
        MINUTES_TO_HOURS = 1 / 60;

        HOURS_TO_SECONDS = 3600;
        SECONDS_TO_HOURS = 1 / 3600;
    end

    methods(Access = public)
        function [obj] = Conversions()
            error("Cannot construct instance of class Conversions.");
        end
    end

end