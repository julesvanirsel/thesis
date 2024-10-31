% user input
direc = fullfile('data', 'paper2');
reload = false; % reload grid
plot_bad_data = true;
save_plot = true;

from_inverions = true;
subtract_background_flow = false;
skip_along_track = false;
do_ephem_fix = true;
is_acc = true;
sfx = strrep(num2str([from_inverions, subtract_background_flow, ...
    skip_along_track, do_ephem_fix, is_acc]), ' ', '');

smooth_distance = 1000e3; % meters
max_v = 4e3; % meters / second
max_dt = 5; % seconds
win = 20; % seconds
outlier_window = 20; % seconds
grd.cell_size = [320, 448];
tdur = 60; % seconds
vvels_cad = 5; % minutes
events = 2;%1:11;

clrs = dictionary(["red", "gre", "blu"], [630, 558, 428]);

%#ok<*UNRCH>

always_reload = false;
if length(events) > 1
    always_reload = true;
end

%% plotting parameters
colorcet = @jules.tools.colorcet;
clm.r = 'L13'; clm.g = 'L14'; clm.b = 'L15'; clm.s = 'L18'; clm.c = 'L19';
scl.v = 1e-3; scl.j = 1e0; scl.c = 1e-3;

% load datetimes, sats, and grid data
evts_fn = fullfile(direc, 'event_data.txt');
evts_num = 0;
fid = fopen(evts_fn);
while ~feof(fid); fgetl(fid); evts_num = evts_num+1; end

i = 1;
times = cell(1, evts_num);
sats = cell(1, evts_num);
arcs = cell(1, evts_num);
sigs = cell(1, evts_num);
vbgs = cell(1, evts_num);
grd.cfgs = cell(1, evts_num);
grd.base.alt_min = 80e3;
grd.base.alt_max = 500e3;
grd.base.alt_scale = [5150 4850 250000 30000];
grd.base.Bincl = 90;
grd.base.lxp = 1;
grd.base.lyp = 1;
fid = fopen(evts_fn);
fgetl(fid);
while ~feof(fid)
    tmp = strsplit(fgetl(fid));
    time = datetime(tmp{2}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
    times{i} = time;
    sats{i} = tmp{3};
    grd.cfgs{i} = grd.base;
    grd.cfgs{i}.times = time;
    grd.cfgs{i}.glat = str2double(tmp{4});
    grd.cfgs{i}.glon = str2double(tmp{5});
    grd.cfgs{i}.xdist = str2double(tmp{6});
    grd.cfgs{i}.ydist = str2double(tmp{7});
    grd.cfgs{i}.x2parms = str2double(strsplit(tmp{8}, ','));
    grd.cfgs{i}.x3parms = str2double(strsplit(tmp{9}, ','));
    arcs{i} = tmp{10};
    sigs{i} = str2double(strsplit(tmp{11}, ','));
    if subtract_background_flow
        vbgs{i} = str2double(strsplit(tmp{12}, ','));
    else
        vbgs{i} = [0,0];
    end
    i = i+1;
end
fclose(fid);

% load skymaps
if not(from_inverions)
    dasc.maph5 = fullfile(direc, 'dasc_data', 'imagery', 'skymap.h5');
    tmp_glat = h5read(dasc.maph5, '/Footpointing/Apex/Latitude');
    tmp_glon = h5read(dasc.maph5, '/Footpointing/Apex/Longitude');
    dasc.glat.red = tmp_glat(1:end-1, :, 3);
    dasc.glat.gre = tmp_glat(1:end-1, :, 2);
    dasc.glat.blu = tmp_glat(1:end-1, :, 1);
    dasc.glat.con = nan(grd.cell_size);
    dasc.glon.red = tmp_glon(1:end-1, :, 3);
    dasc.glon.gre = tmp_glon(1:end-1, :, 2);
    dasc.glon.blu = tmp_glon(1:end-1, :, 1);
    dasc.glon.con = nan(grd.cell_size);
end

% main
for i = events
    % input('Press any key to continue ...')

    time = times{i};
    sat = sats{i};
    vbg = vbgs{i}*scl.v;
    fprintf(pad(' Datetime: %s - Sat: %s ', 80, 'both', '-'), time, sat);
    fprintf('\n\n')

    % load pfisr data
    time.Format = 'uuuuMMdd';
    for pid = 1:99
        try
            pfisr.h5 = dir(fullfile(direc, 'pfisr_data', sprintf('%s.%03i_lp_%imin*lat.h5', time, pid, vvels_cad))).name;
        catch
            continue
        end
        pfisr.h5 = fullfile(direc, 'pfisr_data', pfisr.h5);
        pfisr.time = mean(datetime(h5read(pfisr.h5, '/Time/UnixTime'), 'ConvertFrom', 'posixtime'));
        [~, pfisr.id] = min(abs(pfisr.time - time));
        pfisr.time = pfisr.time(pfisr.id);
        pfisr.dt = pfisr.time - time;
        if abs(seconds(pfisr.dt)) < vvels_cad*60/2
            break
        end
    end
    pfisr.alt = h5read(pfisr.h5, '/VvelsGeoCoords/Altitude');
    [~, pfisr.alt_id] = min(abs(pfisr.alt(1, :) - 110));
    pfisr.glat = h5read(pfisr.h5, '/VvelsGeoCoords/Latitude');
    pfisr.glon = h5read(pfisr.h5, '/VvelsGeoCoords/Longitude') + 360;
    pfisr.glat = pfisr.glat(:, pfisr.alt_id)';
    pfisr.glon = pfisr.glon(:, pfisr.alt_id)';
    pfisr.vvels = h5read(pfisr.h5, '/VvelsGeoCoords/Velocity');
    pfisr.ve = squeeze(pfisr.vvels(1, :, pfisr.alt_id, pfisr.id))*scl.v - vbg(1);
    pfisr.vn = squeeze(pfisr.vvels(2, :, pfisr.alt_id, pfisr.id))*scl.v + vbg(2);
    pfisr.bad = vecnorm([pfisr.ve; pfisr.vn]) > max_v*scl.v; % large flows
    pfisr.ve(pfisr.bad) = nan;
    pfisr.glat_removed = pfisr.glat;
    pfisr.glat_removed(not(pfisr.bad)) = nan;

    % gather dasc data
    if from_inverions
        if is_acc
            dasc.h5 = fullfile(direc, 'dasc_data', sprintf('INV_PKR_%s_%i_ACC.h5', time, second(time, 'secondofday')));
        else
            dasc.h5 = fullfile(direc, 'dasc_data', sprintf('INV_PKR_%s_%i.h5', time, second(time, 'secondofday')));
        end
        dasc.time.con = datetime(h5read(dasc.h5, '/Time/Unix'), 'ConvertFrom', 'posixtime');
        dasc.dt.con = dasc.time.con - time;
        dasc.glat.con = h5read(dasc.h5, '/Coordinates/Latitude');
        dasc.glon.con = h5read(dasc.h5, '/Coordinates/Longitude');
        dasc.SIGP = h5read(dasc.h5, '/Derived/Conductance/Pedersen');
        dasc.SIGH = h5read(dasc.h5, '/Derived/Conductance/Hall');
        tmp_times = datetime(h5read(dasc.h5, '/Optical/Times'), 'ConvertFrom', 'posixtime');
        j = 1;
        for cl = ["Red", "Green", "Blue"]
            c = char(cl); c = string(lower(c(1:3)));
            dasc.time.(c) = tmp_times(j);
            dasc.dt.(c) = dasc.time.(c) - time;
            dasc.(c) = h5read(dasc.h5, ['/Optical/', char(cl)])*scl.c;
            dasc.glat.(c) = dasc.glat.con;
            dasc.glon.(c) = dasc.glon.con;
            j = j + 1;
        end
    else
        imgs_root = fullfile(direc, 'dasc_data', 'imagery', string(time));
        for c = ["red", "gre", "blu"]
            tmp_imgs = {dir(fullfile(imgs_root, ['*', num2str(clrs(c)), '.png'])).name};
            tmp_times = NaT(1, length(tmp_imgs));
            for j = 1:length(tmp_imgs)
                tmp = strsplit(tmp_imgs{j}, '_');
                tmp_times(j) = datetime([tmp{2}, tmp{3}], 'InputFormat', 'uuuuMMddHHmmss');
            end
            [~, dasc.id.(c)] = min(abs(tmp_times - time));
            dasc.(c) = imread(fullfile(imgs_root, tmp_imgs{dasc.id.(c)}));
            dasc.time.(c) = tmp_times(dasc.id.(c));
            dasc.dt.(c) = dasc.time.(c) - time;
            dasc.imgs.(c) = tmp_imgs{dasc.id.(c)};
            dasc.SIGP = nan(grd.cell_size);
            dasc.SIGH = nan(grd.cell_size);
        end
    end

    % gather swarm data
    if skip_along_track
        efi.h5 = sprintf('SW_EXPT_EFI%s*%i%02i%02i*0302_novx.h5', sat, year(time), month(time), day(time));
    else
        efi.h5 = sprintf('SW_EXPT_EFI%s*%i%02i%02i*0302.h5', sat, year(time), month(time), day(time));
    end
    efi.h5 = dir(fullfile(direc, 'swarm_data', efi.h5)).name;
    efi.h5 = fullfile(direc, 'swarm_data', efi.h5);
    efi.glat = h5read(efi.h5, '/GeodeticLatitude110km')';
    efi.glon = wrapTo360(h5read(efi.h5, '/Longitude110km'))';
    efi.times = h5read(efi.h5, '/Timestamp');
    efi.times = datetime(efi.times, 'ConvertFrom', 'posixtime');
    [~, efi.id] = min(abs(efi.times - time));
    efi.ids = (-win:win)*2 + efi.id;
    efi.time = efi.times(efi.id);
    efi.dt = efi.time - time;
    efi.ve = h5read(efi.h5, '/ViE')'*scl.v - vbg(1);
    efi.vn = h5read(efi.h5, '/ViN')'*scl.v - vbg(2); 
    [~, efi.bad] = rmoutliers(vecnorm([efi.ve;efi.vn]), 'movmean', 2*outlier_window); % outliers
    efi.bad = efi.bad | vecnorm([efi.ve;efi.vn]) > max_v*scl.v; % large flows
    if not(plot_bad_data)
        efi.ve(efi.bad) = nan;
    end
    efi.glat_removed = efi.glat;
    efi.glat_removed(not(efi.bad)) = nan;
    
    fac.h5 = sprintf('SW_OPER_FAC%s*%i%02i%02i*0401.h5', sat, year(time), month(time), day(time));
    fac.h5 = dir(fullfile(direc, 'swarm_data', fac.h5)).name;
    fac.h5 = fullfile(direc, 'swarm_data', fac.h5);
    fac.glat = h5read(fac.h5, '/GeodeticLatitude110km')';
    fac.glon = wrapTo360(h5read(fac.h5, '/Longitude110km'))';
    fac.times = h5read(fac.h5, '/Timestamp');
    fac.times = datetime(fac.times, 'ConvertFrom', 'posixtime');
    [~, fac.id] = min(abs(fac.times - efi.time));
    fac.ids = (-win:win) + fac.id;
    fac.time = fac.times(fac.id);
    fac.dt = fac.time - time;
    fac.fac = h5read(fac.h5, '/FAC')*scl.j;
    fac.fac_prec = fac.fac;
    fac.fac_rtrn = fac.fac;
    fac.fac_prec(fac.fac > 0) = nan;
    fac.fac_rtrn(fac.fac <= 0) = nan;
    
    % TII ephemeris fix
    if do_ephem_fix
        fglat = griddedInterpolant(posixtime(fac.times), fac.glat);
        fglon = griddedInterpolant(posixtime(fac.times), fac.glon);
        efi.glat = fglat(posixtime(efi.times));
        efi.glon = fglon(posixtime(efi.times));
    end

    assert(abs(seconds(dasc.dt.red)) < max_dt, 'DASC red time difference: %.3f', seconds(dasc.dt.red))
    assert(abs(seconds(dasc.dt.gre)) < max_dt, 'DASC green time difference: %.3f', seconds(dasc.dt.gre))
    assert(abs(seconds(dasc.dt.blu)) < max_dt, 'DASC blue time difference: %.3f', seconds(dasc.dt.blu))
    assert(abs(seconds(efi.dt)) < max_dt, 'EFI time difference: %.3f', seconds(efi.dt))
    assert(abs(seconds(fac.dt)) < max_dt, 'FAC time difference: %.3f', seconds(fac.dt))
    assert(abs(seconds(pfisr.dt)) < vvels_cad*60/2, 'PFISR time difference: %.3f', seconds(pfisr.dt))
    
    time.Format = 'default';
    fprintf('  Time:                      %s\n', time)
    fprintf('  DASC red:                  %s, dt = %+6.3f seconds\n', dasc.time.red, seconds(dasc.dt.red))
    fprintf('  DASC green:                %s, dt = %+6.3f seconds\n', dasc.time.gre, seconds(dasc.dt.gre))
    fprintf('  DASC blue:                 %s, dt = %+6.3f seconds\n', dasc.time.blu, seconds(dasc.dt.blu))
    fprintf('  EFI: %.3f E, %.3f N @ %s, dt = %+6.3f seconds\n', efi.glon(efi.id), efi.glat(efi.id), efi.time, seconds(efi.dt))
    fprintf('  FAC: %.3f E, %.3f N @ %s, dt = %+6.3f seconds\n', fac.glon(fac.id), fac.glat(fac.id), fac.time, seconds(fac.dt))
    fprintf('  PFISR:                     %s, dt = %+6.3f seconds\n\n', pfisr.time, seconds(pfisr.dt))

    % gather grid data
    grd.cfg = grd.cfgs{i};
    if ~isfield(grd, 'xg') || reload || always_reload
        grd.xg = gemini3d.grid.cartesian(grd.cfg);
        grd.xg = jules.tools.shrink(grd.xg);
        if not(always_reload)
            reload = false;
        end
    else
        warning('Grid is already loaded')
    end

    grd.glat = squeeze(grd.xg.glat(end, :, :));
    grd.glon = squeeze(grd.xg.glon(end, :, :));
    grd.glat = [grd.glat(1, :), grd.glat(:, end)', grd.glat(end, end:-1:1), grd.glat(end:-1:1, 1)'];
    grd.glon = [grd.glon(1, :), grd.glon(:, end)', grd.glon(end, end:-1:1), grd.glon(end:-1:1, 1)'];
    grd.smooth = smooth_distance / range(grd.xg.x2);
    jules.tools.grid_params(grd.cell_size, [grd.cfg.xdist, grd.cfg.ydist]/1e3);
    angle = atan2(grd.xg.glat(1, end, 1) - grd.xg.glat(1, 1, 1), ...
        grd.xg.glon(1, end, 1) - grd.xg.glon(1, 1, 1));

    assert(grd.xg.lx(2) == grd.cell_size(1), 'x2 size not equal to target size.')
    assert(grd.xg.lx(3) == grd.cell_size(2), 'x3 size not equal to target size.')

    % get color limits
    tmp_gre = dasc.gre(inpolygon(dasc.glon.gre(:), dasc.glat.gre(:), grd.glon(:), grd.glat(:)));
    lim.c = [quantile(tmp_gre, 0.05), quantile(tmp_gre, 0.95)];
    tmp_sig = dasc.SIGP(inpolygon(dasc.glon.con(:), dasc.glat.con(:), grd.glon(:), grd.glat(:)));
    lim.p = [quantile(tmp_sig, 0.05), quantile(tmp_sig, 0.95)];
    tmp_sig = dasc.SIGH(inpolygon(dasc.glon.con(:), dasc.glat.con(:), grd.glon(:), grd.glat(:)));
    lim.h = [quantile(tmp_sig, 0.05), quantile(tmp_sig, 0.95)];
    lim.p = lim.p + [-1,1]*range(lim.p)*0.1;
    lim.h = lim.h + [-1,1]*range(lim.h)*0.1;

    if arcs{i} == 'p'
        dasc.arc = smoothdata2(dasc.SIGP, 'gaussian', grd.smooth);
    elseif arcs{i} == 'h'
        dasc.arc = smoothdata2(dasc.SIGH, 'gaussian', grd.smooth);
    end
    dasc.arc(not(inpolygon(dasc.glon.con, dasc.glat.con, grd.glon, grd.glat))) = nan;
    ctrs.A = contour(dasc.glon.con, dasc.glat.con, dasc.arc, [1, 1]*sigs{i}(1));
    ctrs.B = contour(dasc.glon.con, dasc.glat.con, dasc.arc, [1, 1]*sigs{i}(2));
    close all
    try
        [bdry.A.lon, bdry.A.lat] = jules.tools.unpack_contours(ctrs.A, rank=1, min_diff=[2, 0]);
    catch
        bdry.A.lon = nan;
        bdry.A.lat = nan;
    end
    try
        [bdry.B.lon, bdry.B.lat] = jules.tools.unpack_contours(ctrs.B, rank=-1, min_diff=[2, 0]);
    catch
        bdry.B.lon = nan;
        bdry.B.lat = nan;
    end

    lim.x = [min([grd.glon, efi.glon(efi.ids), fac.glon(fac.ids)]), ...
        max([grd.glon, efi.glon(efi.ids), fac.glon(fac.ids)])] + [-1, 1]*0.2;
    lim.y = [min([grd.glat, efi.glat(efi.ids), fac.glat(fac.ids)]), ...
        max([grd.glat, efi.glat(efi.ids), fac.glat(fac.ids)])] + [-1, 1]*0.1;
    % ar = [range(lim.x)*cosd(mean(lim.y)), range(lim.y), 1];
    ar = [1.2, 1, 1];
    
    % plot
    close all
    figure('Position', [10, 60, 1500, 800], 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 10, 10])
    tlo = tiledlayout(2, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
    title(tlo, sprintf([ ...
        '\nTIME: %s UT (UTsec0 = %i)    SAT: %s    GRID: %.0fx%.0fkm    EFI: %+.1fs    FAC: %+.1fs    PFISR: %+.1fs    DLON: %+.2f°    DLAT: %+.2f°    DIR: %s\n' ...
        'BG FLOW: %i,%ikm/s    DASC: %s    PFISR: %s\n' ...
        'EFI: %s    FAC: %s\n'], ...
        time, second(time, 'secondofday')-tdur, sat, range(grd.xg.x2(3:end-2))/1e3, range(grd.xg.x3(3:end-2))/1e3, ...
        seconds(efi.dt), seconds(fac.dt), seconds(pfisr.dt), ...
        fac.glon(fac.id) - efi.glon(efi.id), fac.glat(fac.id) - efi.glat(efi.id), ...
        direc, vbg(1), vbg(2), dasc.h5(length(direc)+1:end), pfisr.h5(length(direc)+1:end), ...
        efi.h5(length(direc)+1:end), fac.h5(length(direc)+1:end)), ...
        'Interpreter', 'none', 'FontSize', 11)

    for j = 1:6
        nexttile
        hold on
        switch j
            case 1
            pcolor(dasc.glon.red, dasc.glat.red, dasc.red)
            quiver(lim.x(1)+0.2, lim.y(1)+0.1, 1, 0, '.-g')
            text(lim.x(1)+0.2, lim.y(1)+0.2, '2 km/s', 'Color', 'g')
            text(lim.x(1)+0.2, lim.y(1)+0.35, 'O removed', 'Color', 'r')
            colormap(gca, colorcet(clm.r))
            clb_label = sprintf('630 nm @ t%+.0f s (kR)', seconds(dasc.dt.red));
            case 2
            pcolor(dasc.glon.gre, dasc.glat.gre, dasc.gre)
            colormap(gca, colorcet(clm.g))
            clb_label = sprintf('558 nm @ t%+.0f s (kR)', seconds(dasc.dt.gre));
            case 3
            pcolor(dasc.glon.blu , dasc.glat.blu , dasc.blu)
            colormap(gca, colorcet(clm.b))
            clb_label = sprintf('428 nm @ t%+.0f s (kR)', seconds(dasc.dt.blu));
            case 4
            pcolor(dasc.glon.con , dasc.glat.con , dasc.SIGP)
            plot(bdry.A.lon, bdry.A.lat, 'k', 'LineWidth', 2)
            plot(bdry.B.lon, bdry.B.lat, '--k', 'LineWidth', 2)
            quiver(fac.glon(fac.ids), fac.glat(fac.ids), ...
                fac.fac_prec(fac.ids)*cos(angle), fac.fac_prec(fac.ids)*sin(angle), 0, '.-b')
            quiver(fac.glon(fac.ids), fac.glat(fac.ids), ...
                fac.fac_rtrn(fac.ids)*cos(angle), fac.fac_rtrn(fac.ids)*sin(angle), 0, '.-r')
            quiver(lim.x(1)+0.2, lim.y(1)+0.1, 1, 0, '.-b')
            text(lim.x(1)+0.2, lim.y(1)+0.2, '1 uA/m^2', 'Color', 'b')
            text(lim.x(1)+0.2, lim.y(1)+0.35, 'Red = Parallel', 'Color', 'r')
            colormap(gca, colorcet(clm.s))
            clim(lim.p)
            clb_label = 'Pedersen Conductance (S)';
            case 5
            pcolor(dasc.glon.con , dasc.glat.con , dasc.SIGH)
            plot(bdry.A.lon, bdry.A.lat, 'k', 'LineWidth', 2)
            plot(bdry.B.lon, bdry.B.lat, '--k', 'LineWidth', 2)
            quiver(fac.glon(fac.ids), fac.glat(fac.ids), ...
                fac.fac_prec(fac.ids)*cos(angle), fac.fac_prec(fac.ids)*sin(angle), 0, '.-b')
            quiver(fac.glon(fac.ids), fac.glat(fac.ids), ...
                fac.fac_rtrn(fac.ids)*cos(angle), fac.fac_rtrn(fac.ids)*sin(angle), 0, '.-r')
            colormap(gca, colorcet(clm.s))
            clim(lim.h)
            clb_label = 'Hall Conductance (S)';
            case 6
            contour(dasc.glon.con , dasc.glat.con , dasc.arc, 6, 'LineWidth', 1)
            plot(bdry.A.lon, bdry.A.lat, 'k', 'LineWidth', 2)
            plot(bdry.B.lon, bdry.B.lat, '--k', 'LineWidth', 2)
            colormap(gca, colorcet('C2'))
            if arcs{i} == 'p'
                clb_label = 'Pedersen Conductance (S)';
            else
                clb_label = 'Hall Conductance (S)';
            end
        end
        if any(j==1:3)
            quiver(efi.glon(efi.ids), efi.glat(efi.ids), efi.ve(efi.ids), efi.vn(efi.ids), 0, '.-g')
            % if skip_along_track
            %     quiver(efi.glon(efi.ids), efi.glat(efi.ids), ...
            %         efi.ve(efi.ids)*cos(angle) - efi.vn(efi.ids)*sin(angle), ...
            %         efi.ve(efi.ids)*sin(angle) + efi.vn(efi.ids)*cos(angle), 0, '.-g')
            % else
            %     quiver(efi.glon(efi.ids), efi.glat(efi.ids), efi.ve(efi.ids), efi.vn(efi.ids), 0, '.-g')
            % end
            scatter(efi.glon(efi.ids), efi.glat_removed(efi.ids), 'or')
            plot(efi.glon(efi.ids), efi.glat(efi.ids), 'w')
            quiver(pfisr.glon, pfisr.glat, pfisr.ve, pfisr.vn, 0, '.-g')
            scatter(pfisr.glon, pfisr.glat_removed, 'or')
            plot(pfisr.glon, pfisr.glat, 'w')
            scatter(efi.glon(efi.id), efi.glat(efi.id), 500, 'xw')
            plot(grd.glon, grd.glat, '--w')
            set(gca, 'Color', [1, 1, 1]/2)
            clim(lim.c)
        else
            plot(fac.glon(fac.ids), fac.glat(fac.ids), 'k')
            scatter(fac.glon(fac.id), fac.glat(fac.id), 500, 'xk')
            plot(grd.glon, grd.glat, '--k')
        end
        shading flat
        clb = colorbar;
        clb.Label.String = clb_label;
        if any(j==4:6); xlabel('geodetic longitude (deg)'); end
        if any(j==[1, 4]); ylabel('geodetic latitude (deg)'); end
        xlim(lim.x); ylim(lim.y)
        pbaspect(ar)
    end
    
    % inset
    outline_glon = dasc.glon.red(not(isnan(dasc.blu)));
    outline_glat = dasc.glat.red(not(isnan(dasc.blu)));
    outline_hull = convhull(outline_glon, outline_glat);
    axes('Position', [-0.015, .875, .09, .09])
    hold on
    plot(outline_glon(outline_hull), outline_glat(outline_hull), 'k')
    plot(grd.glon, grd.glat, 'k')
    plot(efi.glon(efi.ids), efi.glat(efi.ids), 'k')
    plot(pfisr.glon, pfisr.glat, 'k')
    pbaspect([1, 1, 1])
    box off
    set(gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None')
    
    if save_plot
        filename = sprintf('%i%02i%02i_%i_%s_%s.png', ...
            year(time), month(time), day(time), second(time, 'secondofday'), sat, sfx);
        fprintf('Saving %s\n', fullfile(direc, filename))
        exportgraphics(gcf, fullfile(direc, filename), 'Resolution', 600)
        close all
    end
end

fclose all;

%% temp
% t0 = datetime(2023, 2, 10, 9, 50, 0);
% t1 = datetime(2023, 2, 10, 9, 53, 0);
% 
% load data\swop\flow_data_alex\SW_EXPT_EFIA_TCT02_20230210T094600_horiz_geodet.mat
% fac_alex = SW_EXPT_EFIA_TCT02_20230210T094600.FAC_micro;
% glat110_alex = SW_EXPT_EFIA_TCT02_20230210T094600.Footpoint_Lat;
% glon110_alex = SW_EXPT_EFIA_TCT02_20230210T094600.Footpoint_Lon;
% time_alex = NaT(1, length(fac_alex));
% for i = 1:length(fac_alex)
%     tmp = pad(num2str(SW_EXPT_EFIA_TCT02_20230210T094600.UTtime(i)), 10, 'left', '0');
%     time_alex(i) = datetime(tmp, 'InputFormat', 'HHmmss.SSS') - datetime('today') + datetime(2023, 2, 10);
% end
% ids_alex = (time_alex >= t0) & (time_alex <= t1);
% 
% fn_swrm = fullfile('data', 'paper2', 'swarm_data', 'SW_OPER_FACATMS_2F_20230210T000000_20230210T235959_0401.h5');
% fac_swrm = h5read(fn_swrm, '/FAC');
% glat110_swrm = h5read(fn_swrm, '/GeodeticLatitude110km');
% glon110_swrm = h5read(fn_swrm, '/Longitude110km');
% time_swrm = datetime(h5read(fn_swrm, '/Timestamp'), 'ConvertFrom', 'posixtime');
% ids_swrm = (time_swrm >= t0) & (time_swrm <= t1);
% 
% close all
% 
% figure
% hold on
% plot(time_alex(ids_alex)-seconds(1), fac_alex(ids_alex), 'DisplayName', 'Alex (delayed 1 second, or 1 index)')
% plot(time_swrm(ids_swrm), fac_swrm(ids_swrm), 'DisplayName', 'Jules')
% legend
% grid on
% yticks(0)
% xlabel('UT')
% ylabel('Swarm A FAC (uA/m^2)')
% 
% figure
% hold on
% plot(glat110_alex(ids_alex), fac_alex(ids_alex), 'DisplayName', 'Alex')
% plot(glat110_swrm(ids_swrm), fac_swrm(ids_swrm), 'DisplayName', 'Jules')
% legend
% grid on
% yticks(0)
% xlabel('Geodetic Latitude Footpointed to 110 km (deg)')
% ylabel('Swarm A FAC (uA/m^2)')
% 
% figure
% hold on
% plot(glon110_alex(ids_alex), fac_alex(ids_alex), 'DisplayName', 'Alex')
% plot(glon110_swrm(ids_swrm)+0.081, fac_swrm(ids_swrm), 'DisplayName', 'Jules')
% legend
% grid on
% yticks(0)
% xlabel('Geodetic Longitude Footpointed to 110 km (deg)')
% ylabel('Swarm A FAC (uA/m^2)')

%%
% load data\paper2\dasc_data\skymap.mat
% load data\paper2\dasc_data\20230210_35488_dasc.mat
% h5fn = 'data\paper2\dasc_data\skymap.h5';
% h5writeatt(h5fn, '/', 'Site', 'Poker Flat')
% h5writeatt(h5fn, '/', 'Site Altitude', '213 meters')
% h5writeatt(h5fn, '/', 'Site Geodetic Latitude', '213 meters')
% 
% h5create(h5fn, '/Azimuth', size(az), 'Datatype', 'double')
% h5write(h5fn, '/Azimuth', az)
% h5writeatt(h5fn, '/Azimuth', 'Units', 'degrees')
% h5writeatt(h5fn, '/Azimuth', 'Description', 'Azimuth')
% 
% glat = magnetic_footpointing.('110km').lat;
% glon = magnetic_footpointing.('110km').lon;
% sigh = data.green_rayleighs;
% 
% close all
% hold on
% pcolor(glon, glat, sigh)
% quiver(glon110_alex, glat110_alex, fac_alex, 0*fac_alex, '.-g')
% xlim([min(glon(:)), max(glon(:))])
% ylim([min(glat(:)), max(glat(:))])
% shading flat
