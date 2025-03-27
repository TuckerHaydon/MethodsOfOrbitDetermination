function [bool] = mustBeInZeroTo2Pi(x)
    bool = (x >= 0) & (x <= 2*pi);
end