function [x, y] = unpack_contours(cntr,opts)
arguments
    cntr (2,:) double
    opts.rank (1,1) int16 = 1
    opts.min_diff (1,2) double {mustBeNonnegative} = [1,1]
end
%#ok<*AGROW>

mean_ys = [];
lc = size(cntr,2);
i = 1;
j = 1;
while true
    N = cntr(2,i);
    x = cntr(1,i+1:i+N);
    y = cntr(2,i+1:i+N);
    i = i+1+N;
    if i >= lc
        break
    end
    curve = [x; y];
    if all(abs(curve(:,end) - curve(:,1)) < opts.min_diff)
        continue % no closed curves
    end
    mean_ys = [mean_ys,mean(y)];
    curves{j} = curve;
    j = j+1;
end

[~,ids] = sort(mean_ys);
curves_sorted = curves(ids);
if opts.rank == -1
    opts.rank = length(mean_ys);
end
curve = curves_sorted{opts.rank};
x = curve(1,:);
y = curve(2,:);
