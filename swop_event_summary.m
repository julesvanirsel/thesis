% user input
direc = fullfile('data', 'paper2');
reload = true; % reload grid
always_reload = false;
plot_bad_data = true;
save_plot = 1;
events = 1:11;

from_inverions = true;
subtract_background_flow = false;
skip_along_track = false;
do_ephem_fix = true;
is_acc = true;
filter_pfisr = true;
sfx = strrep(num2str([from_inverions, subtract_background_flow, ...
    skip_along_track, do_ephem_fix, is_acc]), ' ', '');

smooth_distance = 5000e3; % meters
max_v = 4e3; % meters / second
max_v_pfisr = 590; % meters / second
max_dt = 5; % seconds
win = 20; % seconds
outlier_window = 20; % seconds
grd.cell_size = [288, 384];
tdur = 60; % seconds
vvels_cad = 5; % minutes

clrs = dictionary(["red", "gre", "blu"], [630, 558, 428]);

%#ok<*UNRCH>

if length(events) > 1
    always_reload = true;
end

% plotting parameters
clm.r = 'L13'; clm.g = 'L14'; clm.b = 'L15'; clm.s = 'L18'; clm.c = 'L19';
scl.v = 1e-3; scl.j = 1e0; scl.c = 1e-3;
fntn = 'Arial';
fnts = 11 * 2;
linw = 1.5;
clbg = [8, 79, 106] / 255;
cltx = [1, 1, 1];
reset(0)
set(0, 'defaultSurfaceEdgeColor', 'flat')
set(0, 'defaultLineLineWidth', linw)
set(0, 'defaultQuiverLineWidth', linw)
jules.tools.setall(0, 'FontName', fntn)
jules.tools.setall(0, 'FontSize', fnts)
jules.tools.setall(0, 'Multiplier', 1)
colorcet = @jules.tools.colorcet;

%% load datetimes, sats, and grid data
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
    pfisr.evvels = h5read(pfisr.h5, '/VvelsGeoCoords/errVelocity');
    pfisr.ve = squeeze(pfisr.vvels(1, :, pfisr.alt_id, pfisr.id))*scl.v - vbg(1);
    pfisr.vn = squeeze(pfisr.vvels(2, :, pfisr.alt_id, pfisr.id))*scl.v + vbg(2);
    pfisr.eve = squeeze(pfisr.evvels(1, :, pfisr.alt_id, pfisr.id))*scl.v;
    pfisr.evn = squeeze(pfisr.evvels(2, :, pfisr.alt_id, pfisr.id))*scl.v;
    pfisr.bad = vecnorm([pfisr.ve; pfisr.vn]) > max_v_pfisr*scl.v; % large flows
    pfisr.ve(pfisr.bad) = nan;
    pfisr.vn(pfisr.bad) = nan;
    pfisr.glat_removed = pfisr.glat;
    pfisr.glat_removed(not(pfisr.bad)) = nan;
    
    %% load superdarn data
    sdarn.h5 = fullfile(direc, 'superdarn_data', 'superdarn_pkr_jules.hdf5');
    sdarn.time = h5read(sdarn.h5, '/df/block0_values');
    sdarn.time = datetime(sdarn.time(:, :)');
    [~, sdarn.id] = min(abs(sdarn.time - time));
    sdarn.time = sdarn.time(sdarn.id);
    sdarn.dt = sdarn.time - time;
    sdarn.b1v = h5read(sdarn.h5, '/df/block1_values');
    sdarn.vme = sdarn.b1v(1, sdarn.id) * scl.v;
    sdarn.vmn = sdarn.b1v(2, sdarn.id) * scl.v;
    sdarn.quality = sdarn.b1v(3, sdarn.id);
    sdarn.nearest_data = sdarn.b1v(4, sdarn.id);

    assert(seconds(sdarn.dt) == 0, 'SuperDARN datetime does not match simulation time.')

    %% gather dasc data
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
    [~, efi.bad] = rmoutliers(vecnorm([efi.ve; efi.vn]), 'movmean', 2*outlier_window); % outliers
    efi.bad = efi.bad | vecnorm([efi.ve; efi.vn]) > max_v*scl.v; % large flows
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
    
    % get background pfisr flow
    [~, ida] = max(grd.glat);
    [~, idb] = max(grd.glon);
    [~, idc] = min(grd.glat);
    [~, idd] = min(grd.glon);
    fpfisr = griddedInterpolant(pfisr.glat, pfisr.glon);
    fbottom = griddedInterpolant([grd.glat(idc), grd.glat(idd)], [grd.glon(idc), grd.glon(idd)]);
    ftop = griddedInterpolant([grd.glat(idb), grd.glat(ida)], [grd.glon(idb), grd.glon(ida)]);
    glat_min = fzero(@(x)(fbottom(x) - fpfisr(x)), 0);
    glat_max = fzero(@(x)(ftop(x) - fpfisr(x)), 0);
    if filter_pfisr
        pfisr.ids = pfisr.glat >= glat_min & pfisr.glat <= glat_max;
        pfisr.ids = pfisr.ids & not(isnan(pfisr.ve));
    else
        pfisr.ids = true(size(pfisr.glat));
    end
    pfisr.glat_avg = pfisr.glat;
    pfisr.glat_avg(not(pfisr.ids)) = nan;
    pfisr.ve_avg = mean(pfisr.ve(pfisr.ids));
    pfisr.vn_avg = mean(pfisr.vn(pfisr.ids));
    mag_dec = atan2(grd.xg.x3(end) - grd.xg.x3(1), grd.xg.x2(end) - grd.xg.x2(1));
    pfisr.vme_avg = pfisr.ve_avg * cos(mag_dec) - pfisr.vn_avg * sin(mag_dec);
    pfisr.vmn_avg = pfisr.ve_avg * sin(mag_dec) + pfisr.vn_avg * cos(mag_dec);

    % rotate sdarn for plotting
    sdarn.ve_avg = sdarn.vmn * sin(mag_dec) + sdarn.vme * cos(mag_dec);
    sdarn.vn_avg = sdarn.vmn * cos(mag_dec) - sdarn.vme * sin(mag_dec);

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
        dasc.arc = smoothdata2(dasc.SIGP, 'gaussian', {grd.smooth, 1});
    elseif arcs{i} == 'h'
        dasc.arc = smoothdata2(dasc.SIGH, 'gaussian', {grd.smooth, 1});
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
    ar = [1, 1, 1];
    
    %% plot
    close all
    figure('Position', [10, 60, 1500, 800], ...
        'PaperUnits', 'inches', 'PaperPosition', [0, 0, 11.5, 7] * 2)
    tlo = tiledlayout(2, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
    title(tlo, sprintf([ ...
        '\nTIME: %s UT (UTsec0 = %i)    SAT: %s    GRID: %.0fx%.0fkm    EFI: %+.1fs    FAC: %+.1fs    PFISR: %+.1fs    DLON: %+.2f°    DLAT: %+.2f°    DIR: %s\n' ...
        'MAG BG FLOW: (%.0f,%.0f) [PFISR] & (%.0f,%.0f) [SDARN] m/s    DASC: %s    PFISR: %s\n' ...
        'EFI: %s    FAC: %s\n'], ...
        time, second(time, 'secondofday')-tdur, sat, range(grd.xg.x2(3:end-2))/1e3, range(grd.xg.x3(3:end-2))/1e3, ...
        seconds(efi.dt), seconds(fac.dt), seconds(pfisr.dt), ...
        fac.glon(fac.id) - efi.glon(efi.id), fac.glat(fac.id) - efi.glat(efi.id), ...
        direc, pfisr.vme_avg / scl.v, pfisr.vmn_avg / scl.v, sdarn.vme / scl.v, sdarn.vmn / scl.v, dasc.h5(length(direc)+1:end), pfisr.h5(length(direc)+1:end), ...
        efi.h5(length(direc)+1:end), fac.h5(length(direc)+1:end)), ...
        'Interpreter', 'none', 'FontSize', fnts * 0.54, 'Color', cltx)

    % legend positions
    legx = lim.x(1) + (lim.x(end) - lim.x(1)) * 0.03;
    legy = lim.y(1) + (lim.y(end) - lim.y(1)) * 0.07;
    legdx = (lim.x(end) - lim.x(1)) / 20;
    legdy = (lim.y(end) - lim.y(1)) / 20;

    for j = 1:6
        ax.(char(64 + j)) = nexttile;
        hold on
        switch j
            case 1
            pcolor(dasc.glon.red, dasc.glat.red, dasc.red)
            quiver(legx, legy - 0.6*legdy, 1, 0, 0, '.-g')
            text(legx, legy + 0*legdy, '1 km/s', 'Color', 'g', 'FontSize', fnts * 0.7)
            text(legx, legy + 1*legdy, 'O removed', 'Color', 'r', 'FontSize', fnts * 0.7)
            text(legx, legy + 2*legdy, 'PFISR BG', 'Color', 'm', 'FontSize', fnts * 0.7)
            text(legx, legy + 3*legdy, 'SDARN BG', 'Color', 'c', 'FontSize', fnts * 0.7)
            % text(quantile(grd.glon, 0.05), quantile(grd.glat, 0.14), ...
                % sprintf('%.0f, %.0f m/s', pfisr.ve_avg * 1e3, pfisr.vn_avg * 1e3), 'Color', 'm', 'FontSize', fnts * 0.7)
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
            quiver(fac.glon(fac.ids), fac.glat(fac.ids), ...
                fac.fac_prec(fac.ids)*cos(angle), fac.fac_prec(fac.ids)*sin(angle), 0, '.-b')
            quiver(fac.glon(fac.ids), fac.glat(fac.ids), ...
                fac.fac_rtrn(fac.ids)*cos(angle), fac.fac_rtrn(fac.ids)*sin(angle), 0, '.-r')
            quiver(lim.x(1)+0.25, lim.y(1) + 0.10, 1, 0, 0, '.-b')
            text(lim.x(1)+0.2, lim.y(1) + 0.22, '1 uA/m^2', 'Color', 'b', 'FontSize', fnts * 0.7)
            text(lim.x(1)+0.2, lim.y(1) + 0.38, 'Red = Parallel', 'Color', 'r', 'FontSize', fnts * 0.7)
            colormap(gca, colorcet(clm.s))
            clim(lim.p)
            clb_label = 'Pedersen Conductance (S)';
            case 5
            pcolor(dasc.glon.con , dasc.glat.con , dasc.SIGH)
            quiver(fac.glon(fac.ids), fac.glat(fac.ids), ...
                fac.fac_prec(fac.ids)*cos(angle), fac.fac_prec(fac.ids)*sin(angle), 0, '.-b')
            quiver(fac.glon(fac.ids), fac.glat(fac.ids), ...
                fac.fac_rtrn(fac.ids)*cos(angle), fac.fac_rtrn(fac.ids)*sin(angle), 0, '.-r')
            colormap(gca, colorcet(clm.s))
            clim(lim.h)
            clb_label = 'Hall Conductance (S)';
            case 6
            contour(dasc.glon.con , dasc.glat.con , dasc.arc, 6, 'LineWidth', 1)
            colormap(gca, colorcet('C2'))
            if arcs{i} == 'p'
                clb_label = 'Pedersen Conductance (S)';
            else
                clb_label = 'Hall Conductance (S)';
            end
        end
        if any(j==1:3)
            quiver(efi.glon(efi.ids), efi.glat(efi.ids), efi.ve(efi.ids), efi.vn(efi.ids), 0, '.-g')
            scatter(efi.glon(efi.ids), efi.glat_removed(efi.ids), 'or', 'LineWidth', linw)
            plot(efi.glon(efi.ids), efi.glat(efi.ids), 'w')
            quiver(pfisr.glon, pfisr.glat, pfisr.ve, pfisr.vn, 0, '.-r')
            quiver(pfisr.glon, pfisr.glat_avg, pfisr.ve, pfisr.vn, 0, '.-g')
            quiver(legx + 2*legdx, legy + 5*legdy, 0, 0, 0 ,'ow')
            quiver(legx + 2*legdx, legy + 5*legdy, pfisr.ve_avg, pfisr.vn_avg, 0 ,'m')
            quiver(legx + 2*legdx, legy + 5*legdy, sdarn.ve_avg, sdarn.vn_avg, 0 ,'c')
            scatter(pfisr.glon, pfisr.glat_removed, 'or', 'LineWidth', linw)
            plot(pfisr.glon, pfisr.glat, 'w')
            scatter(efi.glon(efi.id), efi.glat(efi.id), 500, 'xw', 'LineWidth', linw)
            plot(grd.glon, grd.glat, '--w')
            set(gca, 'Color', [1, 1, 1]/2)
            clim(lim.c)
        else
            plot(fac.glon(fac.ids), fac.glat(fac.ids), 'k')
            scatter(fac.glon(fac.id), fac.glat(fac.id), 500, 'xk', 'LineWidth', linw)
            plot(grd.glon, grd.glat, '--k')
            plot(bdry.A.lon, bdry.A.lat, 'k', 'LineWidth', 2)
            plot(bdry.B.lon, bdry.B.lat, '--k', 'LineWidth', 2)
        end
        shading flat
        clb = colorbar;
        clb.Color = cltx;
        clb.Label.String = clb_label;
        clb.Label.Color = cltx;
        if any(j==4:6); xlabel('geodetic longitude (deg)'); end
        if any(j==[1, 4]); ylabel('geodetic latitude (deg)'); end
        xlim(lim.x); ylim(lim.y)
        pbaspect(ar)
    end
  
    for j = fieldnames(ax)'
        set(ax.(j{1}), 'Color', [1, 1, 1], 'GridColor', cltx, 'MinorGridColor', cltx, ...
            'XColor', cltx, 'YColor', cltx, 'ZColor', cltx)
    end

    % inset
    outline_glon = dasc.glon.red(not(isnan(dasc.blu)));
    outline_glat = dasc.glat.red(not(isnan(dasc.blu)));
    outline_hull = convhull(outline_glon, outline_glat);
    axes('Position', [-0.021, .912, .09, .09])
    hold on
    plot(outline_glon(outline_hull), outline_glat(outline_hull), 'Color', cltx)
    plot(grd.glon, grd.glat, 'Color', cltx)
    plot(efi.glon(efi.ids), efi.glat(efi.ids), 'Color', cltx)
    plot(pfisr.glon, pfisr.glat, 'Color', cltx)
    pbaspect([1, 1, 1])
    box off
    set(gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None')
    
    set(gcf, 'Color', clbg, 'InvertHardcopy', 'off')


    if save_plot
        filename = sprintf('%i%02i%02i_%i_%s_%s.png', ...
            year(time), month(time), day(time), second(time, 'secondofday'), sat, sfx);
        filename = fullfile(direc, filename);
        fprintf('Saving %s\n', filename)
        % exportgraphics(gcf, fullfile(direc, filename), 'Resolution', 600)
        print(gcf, filename, '-dpng', '-r96')
        close all
    end
    
    % plot pfisr data error bars
    figure
    hold on
    plot(pfisr.glat, pfisr.ve, 'r')
    plot(pfisr.glat, pfisr.vn, 'b')
    plot(pfisr.glat, pfisr.ve - pfisr.eve, ':r')
    plot(pfisr.glat, pfisr.ve + pfisr.eve, ':r')
    plot(pfisr.glat, pfisr.vn - pfisr.evn, ':b')
    plot(pfisr.glat, pfisr.vn + pfisr.evn, ':b')
    xlabel('geodetic latitude (°)')
    ylabel('flow (km/s)')
    xlim(lim.y); ylim([-1, 1] * 2 * max_v_pfisr / 1e3)
    legend('geodetic east', 'geodetic north')
    
    if save_plot
        filename = sprintf('%i%02i%02i_%i_%s_%s.png', ...
                year(time), month(time), day(time), second(time, 'secondofday'), sat, sfx);
        filename = fullfile(direc, 'pfisr_data', filename);
        fprintf('Saving %s\n', filename)
        print(gcf, filename, '-dpng', '-r96')
        close all
    end

end

fclose all;
