function [x2, x3] = grid_params(size_cells, size_kms, opts)
arguments
    size_cells (1, 2) int32 {mustBePositive}
    size_kms (1, 2) double {mustBePositive}
    opts.aggression (1, 1) double {mustBeGreaterThan(opts.aggression, 0.1)} = 1
    opts.degdist_factor (1, 1) double {mustBeLessThan(opts.degdist_factor, 1)} = 0.05
    opts.dx_min_factor (1, 1) double {mustBeLessThan(opts.dx_min_factor, 1)} = 0.8
    opts.dx_max_factor (1, 1) double {mustBePositive} = 20
    opts.max_iter (1, 1) int32 {mustBePositive} = 1000
    opts.debug (1, 1) logical = false
end

agg_0 = 1 + (log10(opts.aggression)+1)/10;

xmax_x2 = size_kms(1)/2;
parms_x2(1) = round(xmax_x2*opts.degdist_factor);
parms_x2(2) = round(size_kms(1)/double(size_cells(1))*opts.dx_min_factor, 1);
parms_x2(3) = parms_x2(2)*opts.dx_max_factor;
parms_x2(4) = parms_x2(1);

xmax_x3 = size_kms(2)/2;
parms_x3(1) = round(xmax_x3*opts.degdist_factor);
parms_x3(2) = round(size_kms(2)/double(size_cells(2))*opts.dx_min_factor, 1);
parms_x3(3) = parms_x3(2)*opts.dx_max_factor;
parms_x3(4) = parms_x3(1);

if opts.debug
    close all
    figure
    hold on
    xlabel('Iteration')
    ylabel('Size')
end

if opts.debug
    fprintf('x2_params\nIteration\tLength\tAggression\tx2parms\n')
end
iter = opts.max_iter;
len_old = nan;
agg = agg_0;
parms_x2_0 = parms_x2;
while iter > 0
    x2 = nonuniform_grid(parms_x2, xmax_x2);
    len = length(x2) - 4;
    if opts.debug
        scatter(opts.max_iter-iter, len, 'or')
    end
    if len == len_old
        agg = 1.1*agg;
    end
    if opts.debug
        fprintf('%9i\t%6i\t%10.3f\t%i, %i, %i, %i\n', ...
            opts.max_iter - iter, len, agg, ...
            parms_x2(1)*1e3, parms_x2(2)*1e3, parms_x2(3)*1e3, parms_x2(4)*1e3)
    end
    if len == size_cells(1)
        break
    elseif (len > size_cells(1)) && (parms_x2(1) < 3*parms_x2_0(1))
        parms_x2(1) = round(agg*parms_x2(1), 1);
    elseif len < size_cells(1)
        parms_x2(2) = round(parms_x2(2) / agg, 2);
    else
        parms_x2(1) = parms_x2_0(1);
        parms_x2(4) = round(agg*parms_x2(4), 1);
    end
    len_old = len;
    iter = iter - 1;
end
if iter == 0
    warning('No minimum reached for x2.')
end

if opts.debug
    fprintf('x3_params\nIteration\tLength\tAggression\tx3parms\n')
end
iter = opts.max_iter;
len_old = nan;
agg = agg_0;
parms_x3_0 = parms_x3;
while iter > 0
    x3 = nonuniform_grid(parms_x3, xmax_x3);
    len = length(x3) - 4;
    if opts.debug
        scatter(opts.max_iter-iter, len, 'ob')
    end
    if len == len_old
        agg = 1.1*agg;
    end
    if opts.debug
        fprintf('%9i\t%6i\t%10.3f\t%i, %i, %i, %i\n', ...
            opts.max_iter - iter, len, agg, ...
            parms_x3(1)*1e3, parms_x3(2)*1e3, parms_x3(3)*1e3, parms_x3(4)*1e3)
    end
    if len == size_cells(2)
        break
    elseif (len > size_cells(2)) && (parms_x3(1) < 3*parms_x3_0(1))
        parms_x3(1) = round(agg*parms_x3(1), 1);
    elseif len < size_cells(2)
        parms_x3(2) = round(parms_x3(2) / agg, 2);
    else
        parms_x3(1) = parms_x3_0(1);
        parms_x3(4) = round(agg_0*parms_x3(4), 1);
    end
    len_old = len;
    iter = iter - 1;
end
if iter == 0
    warning('No minimum reached for x3.')
end

if length(x2)-4 ~= size_cells(1)
    error('Target size of x2 is %i, but actual size is %i.', size_cells(1), length(x2)-4)
end
if length(x3)-4 ~= size_cells(2)
    error('Target size of x3 is %i, but actual size is %i.', size_cells(2), length(x3)-4)
end

if opts.debug
    figure
    hold on
    plot(x2(2:end), diff(x2), 'r')
    plot(x3(2:end), diff(x3), 'b')
    xlim([0, inf])
    xlabel('Distance (km)')
    ylabel('Cell size (km)')
    legend('dx2', 'dx3')
end

parms_x2 = parms_x2*1e3;
parms_x3 = parms_x3*1e3;
fprintf('Minimum x2 difference = %i m\n', round(min(diff(x2*1e3))))
fprintf('Maximum x2 difference = %i m\n', round(max(diff(x2*1e3))))
fprintf('Minimum x3 difference = %i m\n', round(min(diff(x3*1e3))))
fprintf('Maximum x3 difference = %i m\n', round(max(diff(x3*1e3))))
fprintf('\nCopy the follong into config.nml:\n\n')
fprintf('x2parms = %i, %i, %i, %i\n', ...
    parms_x2(1), parms_x2(2), parms_x2(3), parms_x2(4))
fprintf('x3parms = %i, %i, %i, %i\n\n', ...
    parms_x3(1), parms_x3(2), parms_x3(3), parms_x3(4))

    function x = nonuniform_grid(parms, xmax)
        degdist = parms(1);    % distance from boundary at which we start to degrade resolution
        dx0 = parms(2);        % min step size for grid
        dxincr = parms(3);     % max step size increase for grid
        ell = parms(4);        % transition length of degradation
        xx2 = xmax-degdist;
        x(1) = 1/2*dx0;
        while x(end) < xmax
            dx = dx0 + dxincr * (1/2+1/2*tanh((x(end)-xx2)/ell));
            x(end+1) = x(end)+dx; %#ok<AGROW>
        end
        x = [-fliplr(x), x];
    end
end