function [bool] = mustBeNonnegativeOrNaN(x)
    bool = (x >= 0) | isnan(x);
end