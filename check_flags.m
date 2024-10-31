direc = fullfile('data', 'paper2');
win = 20; % seconds
scl.v = 1e-3;

event_fn = fullfile(direc, 'event_data.txt');
event_lines = readlines(event_fn)';
num_events = length(event_lines)-2;
for e = 1:num_events
    event_data = strsplit(event_lines{e + 1});
    time = datetime(event_data{2}, 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss''Z');
    sat = event_data{3};
    fprintf('Time: %s    Sat: %s\n', time, sat)
    time.Format = "uuuuMMdd";

    efi.h5 = dir(fullfile(direc, 'swarm_data', sprintf('SW_EXPT_EFI%s*%s*.h5', sat, time))).name;
    efi.h5 = fullfile(direc, 'swarm_data', efi.h5);
    efi.time = datetime(h5read(efi.h5, '/Timestamp'), 'ConvertFrom', 'posixtime');
    efi.cad = 1/mean(rmoutliers(seconds(diff(efi.time))));
    [~, efi.id] = min(abs(efi.time - time));
    efi.ids = efi.id + round(-win*efi.cad:win*efi.cad);
    efi.glat = h5read(efi.h5, '/GeodeticLatitude');
    efi.glon = mod(h5read(efi.h5, '/Longitude'), 360);
    efi.vxh = h5read(efi.h5, '/Vixh')*scl.v;
    efi.vxv = h5read(efi.h5, '/Vixv')*scl.v;
    efi.vy = h5read(efi.h5, '/Viy')*scl.v;
    efi.vz = h5read(efi.h5, '/Viz')*scl.v;
    efi.vxh_err = h5read(efi.h5, '/Vixh_error')*scl.v;
    efi.vxv_err = h5read(efi.h5, '/Vixv_error')*scl.v;
    efi.vy_err = h5read(efi.h5, '/Viy_error')*scl.v;
    efi.vz_err = h5read(efi.h5, '/Viz_error')*scl.v;
    efi.vsatn = h5read(efi.h5, '/VsatN')*scl.v;
    efi.q_flags = int16(h5read(efi.h5, '/Quality_flags'));
    efi.c_flags = int32(h5read(efi.h5, '/Calibration_flags'));

    efi.q_flags = flip(char(pad(string(dec2bin(efi.q_flags)), 16, 'left', '0')), 2);
    efi.c_flags = flip(char(pad(string(dec2bin(efi.c_flags)), 32, 'left', '0')), 2);

    efi.q.vy = str2num(efi.q_flags(:, 3)); %#ok<*ST2NM>
    efi.q.vy(efi.q.vy > 0) = nan;
    efi.c.vxh = str2num(efi.c_flags(:, 1:8));
    efi.c.vxv = str2num(efi.c_flags(:, 9:16));
    efi.c.vy = str2num(efi.c_flags(:, 17:24));
    efi.c.vz = str2num(efi.c_flags(:, 25:32));
    efi.c.vxh(efi.c.vxh == 0) = nan;
    efi.c.vxv(efi.c.vxv == 0) = nan;
    efi.c.vy(efi.c.vy == 0) = nan;
    efi.c.vz(efi.c.vz == 0) = nan;
    efi.c.vxh(efi.c.vxh > 0) = 0;
    efi.c.vxv(efi.c.vxv > 0) = 0;
    efi.c.vy(efi.c.vy > 0) = 0;
    efi.c.vz(efi.c.vz > 0) = 0;

    lps.h5 = dir(fullfile(direc, 'swarm_data', sprintf('SW_OPER_EFI%s*%s*.h5', sat, time))).name;
    lps.h5 = fullfile(direc, 'swarm_data', lps.h5);
    lps.time = datetime(h5read(lps.h5, '/Timestamp'), 'ConvertFrom', 'posixtime');
    lps.cad = 1/mean(rmoutliers(seconds(diff(lps.time))));
    [~, lps.id] = min(abs(lps.time - time));
    lps.ids = lps.id + round(-win*lps.cad:win*lps.cad);
    lps.ne = log10(h5read(lps.h5, '/Ne')*1e6);

    going_north = sign(efi.vsatn(efi.id)) > 0;
    if going_north
        fl = @(x) x;
    else
        fl = @flip;
    end

    close all
    figure('Position', [10, 60, 1520, 820])
    t = tiledlayout(1, 1);
    time.Format = 'default';
    title(t, sprintf('TIME: %s    SAT: %s', time, sat))

    ax1 = axes(t); %#ok<*LAXES>
    ax2 = axes(t);
    for p = 'tl'
        if strcmp(p, 't')
            ax = ax1;
            x = efi.time(efi.ids);
            ax.XTickLabel = string(fl(ax.XTick));
        else
            ax = ax2;
            x = fl(efi.glat(efi.ids));
        end
        hold(ax, 'on')

        plot(ax, x, fl(efi.vxh(efi.ids)), 'DisplayName', 'Vixh', 'Color', [1, 0, 0])
        plot(ax, x, fl(efi.vxh(efi.ids)+efi.vxh_err(efi.ids)), 'DisplayName', 'Vixh+', 'Color', [1, 0, 0], 'LineStyle', ':')
        plot(ax, x, fl(efi.vxh(efi.ids)-efi.vxh_err(efi.ids)), 'DisplayName', 'Vixh-', 'Color', [1, 0, 0], 'LineStyle', ':')
        plot(ax, x, fl(efi.vxv(efi.ids)), 'DisplayName', 'Vixv', 'Color', [1, 1/2, 0])
        plot(ax, x, fl(efi.vxv(efi.ids)+efi.vxv_err(efi.ids)), 'DisplayName', 'Vixv+', 'Color', [1, 1/2, 0], 'LineStyle', ':')
        plot(ax, x, fl(efi.vxv(efi.ids)-efi.vxv_err(efi.ids)), 'DisplayName', 'Vixv-', 'Color', [1, 1/2, 0], 'LineStyle', ':')
        plot(ax, x, fl(efi.vy(efi.ids)), 'DisplayName', 'Viy', 'Color', [0, 0, 1])
        plot(ax, x, fl(efi.vy(efi.ids)+efi.vy_err(efi.ids)), 'DisplayName', 'Viy+', 'Color', [0, 0, 1], 'LineStyle', ':')
        plot(ax, x, fl(efi.vy(efi.ids)-efi.vy_err(efi.ids)), 'DisplayName', 'Viy-', 'Color', [0, 0, 1], 'LineStyle', ':')
        plot(ax, x, fl(efi.vz(efi.ids)), 'DisplayName', 'Viz', 'Color', [0, 1/2, 1/2])
        plot(ax, x, fl(efi.vz(efi.ids)+efi.vz_err(efi.ids)), 'DisplayName', 'Viz+', 'Color', [0, 1/2, 1/2], 'LineStyle', ':')
        plot(ax, x, fl(efi.vz(efi.ids)-efi.vz_err(efi.ids)), 'DisplayName', 'Viz-', 'Color', [0, 1/2, 1/2], 'LineStyle', ':')
        scatter(ax, x, fl(efi.q.vy(efi.ids)), 200, 'DisplayName', 'Bad Qual. Viy', 'MarkerEdgeColor', [1, 0 , 0])
        scatter(ax, x, fl(efi.c.vxh(efi.ids)), 100, 'DisplayName', 'No Cal. Vixh', 'MarkerEdgeColor', [0, 1/2 , 0], 'Marker', 'x')
        scatter(ax, x, fl(efi.c.vxh(efi.ids)), 100, 'DisplayName', 'No Cal. Vixv', 'MarkerEdgeColor', [0, 1/2 , 0], 'Marker', '+')
        scatter(ax, x, fl(efi.c.vxh(efi.ids)), 100, 'DisplayName', 'No Cal. Viy', 'MarkerEdgeColor', [0, 1/2 , 0], 'Marker', '*')
        scatter(ax, x, fl(efi.c.vxh(efi.ids)), 100, 'DisplayName', 'No Cal. Viz', 'MarkerEdgeColor', [0, 1/2 , 0], 'Marker', 'square')
        
        ax.Box = 'off';
        ax.XTick = x(round(length(x)*[0.125, 0.375, 0.625, 0.875]));
        grid(ax)
        xlim(ax, [x(1), x(end)])
        ylim(ax, [-1, 1]*3)
    end

    plot(ax1, lps.time(lps.ids), 2*(lps.ne(lps.ids) - 10), '--k', 'DisplayName', 'Ne')

    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'right';
    ax2.Color = 'none';
    legend(ax1, 'Location', 'northwestoutside')
    % legend(ax2)
    xlabel(ax2, 'Geodetic latitude')
    ylabel(ax1, 'Ion flow (km/s)')
    ylabel(ax2, 'log_{10} Electron Density (m^{-3})')
    ax2.YTickLabel = {"", "9", "", "10", "", "11", ""};

    time.Format = "uuuuMMdd";
    exportgraphics(gcf, fullfile(direc, 'swarm_data', sprintf('%s_%i_%s.png', time, second(time, 'secondofday'), sat)), "Resolution", 600)

    % input('Press any key to continue ... \n')
end

%%
% for f = h5fns
%     h5fn = fullfile(direc, 'swarm_data', f{1});
%     time = datetime(h5read(h5fn, '/Timestamp'), 'ConvertFrom', 'posixtime');
%     cad = 1/seconds(median(diff(time)));
%     glat = h5read(h5fn, '/GeodeticLatitude');
%     glon = mod(h5read(h5fn, '/Longitude'), 360);
%     vxh = h5read(h5fn, '/Vixh')*scl.v;
%     vxv = h5read(h5fn, '/Vixv')*scl.v;
%     vy = h5read(h5fn, '/Viy')*scl.v;
%     vz = h5read(h5fn, '/Viz')*scl.v;
%     vxh_err = h5read(h5fn, '/Vixh_error')*scl.v;
%     vxv_err = h5read(h5fn, '/Vixv_error')*scl.v;
%     vy_err = h5read(h5fn, '/Viy_error')*scl.v;
%     vz_err = h5read(h5fn, '/Viz_error')*scl.v;
%     vsatn = h5read(h5fn, '/VsatN')*scl.v;
%     q_flags = h5read(h5fn, '/Quality_flags');
%     c_flags = h5read(h5fn, '/Calibration_flags');
% 
%     for i = 1:length(event_times)
%         [dt, id] = min(abs(seconds(time - event_times(i))));
%         if dt < 1
%             disp(event_times(i))
%             break
%         end
%     end
%     ids = id + (-win*cad:win*cad);
% 
%     dir = sign(median(vsatn(id)));
%     if dir > 0
%         fl = @(x) x;
%     else
%         fl = @flip;
%     end
%     disp(dir)
% 
%     close all
%     figure('Position', [10, 60, 1520, 820])
% 
%     t = tiledlayout(1, 1);
%     ax1 = axes(t); %#ok<*LAXES>
%     ax2 = axes(t);
%     for p = 'tl'
%         if strcmp(p, 't')
%             ax = ax1;
%             x = time(ids);
%         else
%             ax = ax2;
%             x = fl(glat(ids));
%         end
%         hold(ax, 'on')
% 
%         plot(ax, x, fl(vxh(ids)), 'DisplayName', 'Vixh', 'Color', [1, 0, 0])
%         plot(ax, x, fl(vxv(ids)), 'DisplayName', 'Vixv', 'Color', [1, 1/2, 0])
%         plot(ax, x, fl(vy(ids)), 'DisplayName', 'Viy', 'Color', [0, 0, 1])
%         plot(ax, x, fl(vz(ids)), 'DisplayName', 'Viz', 'Color', [0, 1/2, 1/2])
% 
%         ax.Box = 'off';
%         ax.XTick = x(round(length(x)*[0.125, 0.375, 0.625, 0.875]));
%         if dir < 0 & strcmp(p, 't')
%             ax.XTickLabel = string(ax.XTick);
%         elseif dir > 0 & strcmp(p, 't')
%             ax.XTickLabel = string(flip(ax.XTick));
%         end
%         grid(ax)
%         xlim(ax, [x(1), x(end)])
%         ylim(ax, [-1, 1]*5)
%     end
% 
%     ax2.XAxisLocation = 'top';
%     ax2.Color = 'none';
%     yticks(ax2, [])
%     legend()
% 
%     input('Press any key to continue ... ')
% end

%%


%%
% h5fn = 'data\swop\swarm_data\SW_EXPT_EFIB_TCT02_20230319T040252_20230319T163507_0302.h5';
% mlat = flip(h5read(h5fn, '/MagneticLatitude'));
% ids = h5read(h5fn, '/ids');
% vxh = h5read(h5fn, '/Vixh');
% vxv = h5read(h5fn, '/Vixv');
% vy = h5read(h5fn, '/Viy');
% vz = h5read(h5fn, '/Viz');
% vxhe = h5read(h5fn, '/Vixh_error');
% vxve = h5read(h5fn, '/Vixv_error');
% vye = h5read(h5fn, '/Viy_error');
% vze = h5read(h5fn, '/Viz_error');
% qf = h5read(h5fn, '/Quality_flags');
% cf = h5read(h5fn, '/Calibration_flags');
% vxh = flip(vxh(ids(1):ids(2)-1));
% vxv = flip(vxv(ids(1):ids(2)-1));
% vy = flip(vy(ids(1):ids(2)-1));
% vz = flip(vz(ids(1):ids(2)-1));
% vxhe = flip(vxhe(ids(1):ids(2)-1));
% vxve = flip(vxve(ids(1):ids(2)-1));
% vye = flip(vye(ids(1):ids(2)-1));
% vze = flip(vze(ids(1):ids(2)-1));
% qf = flip(qf(ids(1):ids(2)-1));
% cf = flip(cf(ids(1):ids(2)-1));
% 
% mlat0 = 64.7382;
% mlat1 = 69.8633;
% 
% [~, id0] = min(abs(mlat-mlat0));
% [~, id1] = min(abs(mlat-mlat1));
% 
% close all
% figure(1)
% hold on
% plot(mlat, vxh, 'r')
% plot(mlat, vxhe*100, ':r')
% plot(mlat, vxv, 'b')
% plot(mlat, vxve*100, ':b')
% plot(in_situ.pos(:, 2), in_situ.flow(:, 1))
% plot(in_situ.pos(:, 2), in_situ.flow(:, 2))
% legend('vxh', '100 x vxh_{error}', 'vxv', '100 x vxv_{error}')
% xlim([mlat0, mlat1])
% ylim([-1, 1]*3000)
% 
% figure(2)
% hold on
% plot(mlat, vy, 'r')
% plot(mlat, vye*100, ':r')
% plot(mlat, vz, 'b')
% plot(mlat, vze*100, ':b')
% plot(in_situ.pos(:, 2), in_situ.flow(:, 1))
% plot(in_situ.pos(:, 2), in_situ.flow(:, 2))
% legend('vy', '100 x vy_{error}', 'vz', '100 x vze_{error}')
% xlim([mlat0, mlat1])
% ylim([-1, 1]*3000)


