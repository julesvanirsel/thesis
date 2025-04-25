% direc_root = fullfile('..', '..', 'public_html', 'Gemini3D');
% direc = fullfile(direc_root, 'swop_20230210_35487_AC_09_SD');
% cfg = gemini3d.read.config(direc);
% time = cfg.times(end);

utcs = [[2, 10, 9, 51, 27]; ...
       [2, 12, 10, 22, 11]; ...
       [3, 4, 7, 30, 12]; ...
       [3, 4, 10, 13, 49]; ...
       [3, 14, 6, 49, 7]; ...
       [3, 19, 8, 23, 30]]';

utcs2 = [[2, 10, 9, 50, 00]; ...
       [2, 12, 10, 22, 00]; ...
       [3, 4, 7, 30, 00]; ...
       [3, 4, 10, 12, 00]; ...
       [3, 14, 6, 48, 00]; ...
       [3, 19, 8, 22, 00]]';

close all
fig = figure;
tiledlayout(6, 2, 'TileSpacing', 'tight', 'Padding', 'tight')
set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 6.5, 7]*2, ...
    'Position', [100, 100, 1200, 1200*7/6.5])

for utc = utcs
for it = ["m1s", "m1m"]
    time = datetime([2023, utc']);
    time = time - minutes(10); % 10 minute delay
    % fprintf('\n\n%s\n', time)
    time.Format = 'uuuuMMdd';

    dscvr_direc = fullfile('vanirsel_etal_2025', 'data', 'superdarn');
    dscvr_name = dir(sprintf('%s%s*%s*_s%s000000*.nc', ...
        dscvr_direc, filesep, it, time)).name;
    dscvr_file = fullfile(dscvr_direc, dscvr_name);
    dscvr_time = datetime(1970, 1, 1, 0, 0, 0, ncread(dscvr_file, 'time'));
    cad = diff(posixtime(dscvr_time));
    idl = round(15*60/cad);
    [~, id] = min(abs(time - dscvr_time));
    ids = id + (-idl:idl);

    bt = ncread(dscvr_file, 'bt'); % total
    bx = ncread(dscvr_file, 'bx_gsm'); % Geocentric Solar Ecliptic / Magnetic
    by = ncread(dscvr_file, 'by_gsm');
    bz = ncread(dscvr_file, 'bz_gsm');

    angl = rad2deg(atan2(bz, by));

    % bt_gse = vecnorm([dscvr_bx_gse, dscvr_by_gse, dscvr_bz_gse]');

    % fprintf('Bt = %.2f nT\n\n', dscvr_bt(id))
    % 
    % fprintf('GSE Bx = %.2f nT\n', dscvr_bx_gse(id))
    % fprintf('GSE By = %.2f nT\n', dscvr_by_gse(id))
    % fprintf('GSE Bz = %.2f nT\n', dscvr_bz_gse(id))

    % fprintf('GSM Bx = %.2f nT\n', dscvr_bx_gsm(id))
    % fprintf('GSM By = %.2f nT\n', dscvr_by_gsm(id))
    % fprintf('GSM Bz = %.2f nT\n', dscvr_bz_gsm(id))
    nexttile
    hold on
    time.Format = 'uuuu-MM-dd'' ''HH:mm:ss';
    title(sprintf('%s (GSM %s)', time, it))

    yyaxis left
    plot(dscvr_time(ids), bt(ids), 'k')
    plot(dscvr_time(ids), bx(ids))
    plot(dscvr_time(ids), by(ids))
    plot(dscvr_time(ids), bz(ids))
    ylim([-6, 14])

    yyaxis right
    plot(dscvr_time(ids), angl(ids))
    legend( ...
        sprintf('t (%.1f nT', median(bt(ids))), ...
        sprintf('x (%.1f nT)', median(bx(ids))), ...
        sprintf('y (%.1f nT)', median(by(ids))), ...
        sprintf('z (%.1f nT)', median(bz(ids))), ...
        sprintf('a (%.0f°)', median(angl(ids))), ...
        'location', 'eastoutside')
    xlim([dscvr_time(ids(1)), dscvr_time(ids(end))])
    % ylim([-6, 14])
    yticks(-180:180:180)
    ylim([-180, 180])
    xticks([time - minutes(10), time, time + minutes(10)])
    grid on
end
end

print(fig, 'swop_imf.png', '-dpng', '-r96');
close all
