function [bool] = mustBeInZeroToPi(x)
    bool = (x >= 0) & (x <= pi);
end