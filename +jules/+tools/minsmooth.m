function [y,w] = minsmooth(x)
    w = 0;
    is_unique = false;
    while not(is_unique)
        w = w+1;
        y = smooth(x,w);
        is_unique = numel(y) == numel(unique(y));
    end
end