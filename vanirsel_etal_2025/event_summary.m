direc_root = fullfile('..', '..', '..', 'public_html', 'Gemini3D');
direcs = string(direc_root) + filesep + [ ...
    "swop_20230210_35487_AC_09_SD"; ...
    "swop_20230212_37331_C_09_SD"; ...
    "swop_20230304_27012_C_09_SD"; ...
    "swop_20230304_36829_B_09_SD"; ...
    "swop_20230314_24547_AC_09_SD"; ...
    "swop_20230319_30210_B_09_SD"; ...
    ]';

scl.x = 1e-3;   unt.x = 'km';
scl.c = 1e-3;   unt.c = 'keV';    clm.c = 'L17';
scl.j = 1e6;    unt.j = 'uA/m^2'; clm.j = 'D1A';
scl.U = 1e+3;   unt.U = 'mW/m^2'; clm.U = 'L19';
scl.jq = scl.j / 2;

reset(0)
setall(0, 'FontName', 'Arial')
setall(0, 'FontSize', 10*2)
setall(0, 'Multiplier', 1)
set(0, 'defaultAxesFontSizeMode','manual')
set(0, 'defaultSurfaceEdgeColor','flat')

close all
fig = figure;
tlo = tiledlayout(6, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 6.5, 7]*2, ...
    'Position', [100, 100, 1200, 1200*7/6.5])

for i = 1:numel(direcs)
    direc = direcs(i);
    cfg = gemini3d.read.config(direc);
    time = cfg.times(end);
    filename = gemini3d.datelab(time) + '.h5';
    time.Format = 'MMM d, H:mm';
    disp(time)
    
    mlon = h5read(fullfile(direc, cfg.prec_dir, 'simgrid.h5'), '/mlon');
    mlat = h5read(fullfile(direc, cfg.prec_dir, 'simgrid.h5'), '/mlat');
    [MLON, MLAT] = ndgrid(mlon, mlat);

    x2 = h5read(fullfile(direc, 'inputs', 'simgrid.h5'), '/x2');
    x3 = h5read(fullfile(direc, 'inputs', 'simgrid.h5'), '/x3');
    bdry.A = h5read(fullfile(direc, 'ext', 'current.h5'), '/Boundary/Primary');
    bdry.B = h5read(fullfile(direc, 'ext', 'current.h5'), '/Boundary/Secondary');
    bdry.A(1, :) = interp1(x2(3:end-2), mlon, bdry.A(1, :));
    bdry.A(2, :) = interp1(x3(3:end-2), mlat, bdry.A(2, :));
    bdry.B(1, :) = interp1(x2(3:end-2), mlon, bdry.B(1, :));
    bdry.B(2, :) = interp1(x3(3:end-2), mlat, bdry.B(2, :));

    sats = [h5info(fullfile(direc, 'ext', 'tracks.h5')).Groups.Name];
    sats = strrep(sats, '/', '');
    tracks = struct;
    for sat = sats
        tracks.(sat).x2 = h5read(fullfile(direc, 'ext', 'tracks.h5'), ...
            sprintf('/%s/Coordinates/Magnetic/East', sat));
        tracks.(sat).x3 = h5read(fullfile(direc, 'ext', 'tracks.h5'), ...
            sprintf('/%s/Coordinates/Magnetic/North', sat));
        tracks.(sat).fac = h5read(fullfile(direc, 'ext', 'tracks.h5'), ...
            sprintf('/%s/Current/FieldAligned', sat))' * scl.jq;
        tracks.(sat).mlon = interp1(x2(3:end-2), mlon, tracks.(sat).x2);
        tracks.(sat).mlat = interp1(x3(3:end-2), mlat, tracks.(sat).x3);
        % tmp = tracks.(sat).fac(tracks.(sat).x3 <= max(x3) & tracks.(sat).x3 >= min(x3));
        % fprintf('FAC ranges from %.1f to %.1f for swarm %s on %s\n', min(tmp), max(tmp), sat, time)
    end

    Qp = h5read(fullfile(direc, cfg.prec_dir, filename), '/Qp'); % already has mW units
    Ud = h5read(fullfile(direc, cfg.prec_dir, filename), '/E0p') * scl.c;
    j1 = h5read(fullfile(direc, cfg.E0_dir, filename), '/Vmaxx1it') * scl.j;
    % fprintf('FAC ranges from %.1f to %.1f on %s\n', quantile(-j1(:), 0.1), quantile(-j1(:), 0.9), time)

    lim.x = [min(mlon), max(mlon)];
    lim.y = [min(mlat), max(mlat)];
    lim.j = [-1, 1] * quantile(abs(j1(:)), 0.99);

    nexttile
    hold on
    pcolor(MLON, MLAT, Qp)
    plot(bdry.A(1, :), bdry.A(2, :), 'k')
    plot(bdry.B(1, :), bdry.B(2, :), '--k')
    for sat = sats
        quiver(tracks.(sat).mlon, tracks.(sat).mlat, ...
            tracks.(sat).fac, zeros(size(tracks.(sat).fac)), 0, '.-k')
    end
    colormap(gca, colorcet(clm.U))
    colorbar
    xlim(lim.x); ylim(lim.y)
    ylabel('Mag. Lat. (°)')
    ytickformat(gca, '%4.1f')
    if i == 1
        title('Energy Flux, Q_p (mW/m^2)')
    elseif i == 6
        xticks([256, 260, 264])
    end

    nexttile
    hold on
    pcolor(MLON, MLAT, Ud)
    plot(bdry.A(1, :), bdry.A(2, :), 'k')
    plot(bdry.B(1, :), bdry.B(2, :), '--k')
    colormap(gca, colorcet(clm.c))
    colorbar
    xlim(lim.x); ylim(lim.y)
    yticks([])
    if i == 1
        title('Acc. Potential, U_d (keV)')
    elseif i == 6
        xticks([256, 260, 264])
    end

    nexttile
    hold on
    pcolor(MLON, MLAT, -j1)
    plot(bdry.A(1, :), bdry.A(2, :), 'k')
    plot(bdry.B(1, :), bdry.B(2, :), '--k')
    colormap(gca, colorcet(clm.j))
    clb = colorbar;
    clb.Label.String = sprintf('%s UT', time);
    xlim(lim.x); ylim(lim.y); clim(lim.j)
    yticks([])
    if i==1
        title('FAC, j_{||} (\muA/m^2)')
    elseif i == 6
        xticks([256, 260, 264])
    end
end

print(fig, 'plots/00_event_summary.png', '-dpng', '-r96');
close all

ftsz = 11;
ll = 60;
tt = 30;
ww = 378;
hh = 215;
im = imread('plots/00_event_summary.png');
for i = 1:18
    px = ll + ww * mod(i - 1, 3);
    py = tt + hh * floor((i - 1) / 3);
    im = insertText(im, [px, py], char(64 + i), ...
        'AnchorPoint', 'LeftTop', 'BoxOpacity', 0, 'FontSize', ftsz*2);
end
imwrite(im, 'plots/00_event_summary.png', 'png')