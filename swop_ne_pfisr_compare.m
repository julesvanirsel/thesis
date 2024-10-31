direc = fullfile('..', '..', 'public_html', 'Gemini3D', 'swop_0314_AC_02');
int_time = 20; % minutes

pfisr_root = fullfile('data', 'paper2', 'pfisr_data');

cfg = gemini3d.read.config(direc);
xg = gemini3d.read.grid(direc);
xg = jules.tools.shrink(xg);
gem.time = cfg.times(end);
dat = gemini3d.read.frame(direc, 'time', gem.time, 'var', 'ne');

%%
gem.ne = double(dat.ne);
gem.mlat = double(xg.mlat)';
gem.mlon = double(xg.mlon);
gem.alt = double(xg.alt(:,1,1))' / 1e3;
[gem.ALT, gem.MLON, gem.MLAT] = ndgrid(gem.alt, gem.mlon, gem.mlat);
gem.fne = griddedInterpolant(gem.ALT, gem.MLON, gem.MLAT, gem.ne, 'linear', 'none');

%%
gem.time.Format = 'uuuuMMdd';
for c = ["ac", "lp"]
    h5fn = sprintf('%s.*_%s_%imin-fitcal.h5', gem.time, c, int_time);
    h5fn = fullfile(pfisr_root, h5fn);
    h5fn = fullfile(pfisr_root, dir(h5fn).name);

    time = datetime(h5read(h5fn, '/Time/UnixTime'), 'ConvertFrom', 'posixtime');
    time = mean(time);
    [~, pfisr.(c).record_id] = min(abs(time - gem.time));
    pfisr.(c).time = time(pfisr.(c).record_id);

    ne = h5read(h5fn, '/FittedParams/Ne');
    pfisr.(c).ne = double(ne(:, :, pfisr.(c).record_id));
    pfisr.(c).lat = double(h5read(h5fn, '/Geomag/Latitude'));
    pfisr.(c).lon = double(h5read(h5fn, '/Geomag/Longitude'));
    pfisr.(c).alt = double(h5read(h5fn, '/Geomag/Altitude')) / 1e3;
    [pfisr.(c).theta, pfisr.(c).phi] = gemini3d.geog2geomag(pfisr.(c).lat, pfisr.(c).lon);
    pfisr.(c).mlat = 90 - rad2deg(pfisr.(c).theta);
    pfisr.(c).mlon = rad2deg(pfisr.(c).phi);
    num_beams = size(ne, 2);

    pfisr.(c).alt_reg = linspace(min(pfisr.(c).alt(:)), max(pfisr.(c).alt(:)), 128);
    pfisr.(c).int_ne = nan(length(pfisr.(c).alt_reg), num_beams);
    for b = 1:num_beams
        ids = not(isnan(pfisr.(c).alt(:, b)));
        fne.(c) = griddedInterpolant(pfisr.(c).alt(ids, b), pfisr.(c).ne(ids, b), 'linear', 'none');
        pfisr.(c).int_ne(:, b) = fne.(c)(pfisr.(c).alt_reg);
    end
    pfisr.(c).mean_ne_reg = nanmedian(pfisr.(c).int_ne, 2)';

    pfisr.(c).loess_ne_reg = smooth(repmat(pfisr.(c).alt_reg, [1, num_beams]), pfisr.(c).int_ne(:), 100, 'rlowess');
    pfisr.(c).loess_ne_reg = pfisr.(c).loess_ne_reg(1:length(pfisr.(c).alt_reg))';

    gem.(c).int_ne = gem.fne(pfisr.(c).alt, pfisr.(c).mlon, pfisr.(c).mlat);
    gem.(c).mean_ne = nanmedian(gem.(c).int_ne,2)';
end

%%
close all
figure('PaperPosition',[0,0,10,10], 'PaperUnits', 'inches')
hold on
for c = ["ac", "lp"]
    clr = 'b';
    if strcmp(c,"ac")
        clr = 'r';
    end
    plot(log10(pfisr.(c).loess_ne_reg), pfisr.(c).alt_reg, clr, 'LineWidth', 3, 'DisplayName', sprintf('LOESS PFISR (%s)', c))
    scatter(log10(pfisr.(c).ne(:)), pfisr.(c).alt(:), clr, 'Marker', '.', 'DisplayName', sprintf('PFISR (%s)', c))
    scatter(log10(gem.(c).int_ne(:)), pfisr.(c).alt(:) , clr, 'Marker', 'o', 'DisplayName', sprintf('GEMINI (%s)', c))
end
for i = [1, xg.lx(3)/2, xg.lx(3)]
    clr = [0, double(i)/double(xg.lx(3)), 0];
    plot(log10(mean(squeeze(gem.ne(:,:,i)), 2)), gem.alt, 'Color', clr, 'LineWidth', 3, 'DisplayName', sprintf('GEMINI @ %.1f° mag. lat.', gem.mlat(i)))
end
xlabel('log_{10} Electron density (m^{-3})')
ylabel('Altitude (km)')
xlim([7,14])
ylim([80,200])
% legend('PFISR LOESS (AC)', 'PFISR (AC)', 'GEMINI (AC)', ...
    % 'PFISR LOESS (LP)', 'PFISR (LP)', 'GEMINI (LP)')
legend
title(sprintf('%s UT', pfisr.(c).time))
saveas(gcf, 'test.png')