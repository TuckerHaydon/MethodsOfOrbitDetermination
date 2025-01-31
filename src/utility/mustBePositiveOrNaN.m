function [bool] = mustBePositiveOrNaN(x)
    bool = (x > 0) | isnan(x);
end