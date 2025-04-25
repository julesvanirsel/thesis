clbg = [200, 200, 200] / 255;
clr1 = [244, 187, 129] / 255;
clr2 = [100, 255, 255] / 255;
clr3 = [0, 0, 0];
clr4 = [255, 30, 30] / 255;
ftsz = 15;
letter_pos = [[5, 5]; [590, 5]; [5, 640]; [590, 640]];

linw = 2;
reset(0)
setall(0, 'FontName', 'Arial')
setall(0, 'FontSize', 10*2)
setall(0, 'Multiplier', 1)
setall(0, 'LineWidth', 1.5)
set(0, 'defaultAxesFontSizeMode', 'manual')
set(0, 'defaultSurfaceEdgeColor', 'flat')
set(0, 'defaultAxesLineWidth', 1)

clm.U = 'L19'; clm.c = 'L17';

event_fn = fullfile('data', 'event_data.txt');
lines = readlines(event_fn)';
data = strsplit(lines(3));

%%
Q_filter_quantile = 0.6;
red_filter_quantile = 0.3;
time = datetime(data(2), 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss''Z');
time.Format = 'uuuuMMdd';
h5fn = sprintf('INV_PKR_%s_%i.h5', time, second(time, 'secondofday'));
h5fn = fullfile('data', 'dasc', h5fn);
glat = h5read(h5fn, '/Coordinates/Latitude');
glon = h5read(h5fn, '/Coordinates/Longitude');
Q = h5read(h5fn, '/Derived/Energy/Flux');
E0 = h5read(h5fn, '/Derived/Energy/Characteristic') / 1e3;
red = h5read(h5fn, '/Optical/Red');
glat(isnan(E0)) = nan;
glon(isnan(E0)) = nan;

% filter of large Q and dim red
E0_filtered = E0;
E0_filtered(Q > quantile(Q(:), Q_filter_quantile)) = nan;
E0_filtered(red < quantile(red(:), red_filter_quantile)) = nan;

% determine energy histograms
bin_edges = linspace(0.010, 3, 301);
bin_means = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
hist1 = histcounts(E0, bin_edges);
hist2 = histcounts(E0_filtered, bin_edges);

% determine thermal energies by maximum
[m1, id1] = max(hist1);
[m2, id2] = max(hist2);
temp_max_unfiltered = bin_means(id1);
temp_max_filtered = bin_means(id2);

% determine thermal energies through gaussian fits
start1 = [exp(1)*m1 / temp_max_unfiltered^2, temp_max_unfiltered];
start2 = [exp(1)*m2 / temp_max_filtered^2, temp_max_filtered];
[f1, ~] = fit(bin_means', hist1', 'a*x^2*exp(-(x/b)^2)', 'start', start1);
[f2, ~] = fit(bin_means', hist2', 'a*x^2*exp(-(x/b)^2)', 'start', start2);
temp_fit_unfiltered = f1.b;
temp_fit_filtered = f2.b;

%%
close all
fig = figure('PaperUnits', 'inches', 'PaperPosition', [0, 0, 6.5, 6.5]*2, ...
    'Position', [100, 100, 1200, 1200]);
tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding','tight');

ar = [1, 1, 1];

ax1 = nexttile;
pcolor(glon, glat, Q)
shading flat
colorbar;
colormap(gca, colorcet(clm.U))
clim([0, 21])
xlabel('Geodetic longitude (°)')
ylabel('Geodetic latitude (°)')
pbaspect(ar)
title('Total energy flux, Q_p (mW/m^2)')

ax2 = nexttile;
pcolor(glon, glat, E0)
shading flat
colorbar;
colormap(gca, colorcet(clm.c))
clim([0, 3.1])
xlabel('Geodetic longitude (°)')
ylabel('Geodetic latitude (°)')
pbaspect(ar)
title('Unaccelerated characteristic energy, E_0 (keV)')

ax3 = nexttile;
pcolor(glon, glat, E0_filtered)
shading flat
clb = colorbar;
colormap(gca, colorcet(clm.c))
clim([0, 3.1])
xlabel('Geodetic longitude (°)')
ylabel('Geodetic latitude (°)')
pbaspect(ar)
title('Filtered E_0 (keV)')

nexttile;
hold on
histogram(E0, bin_edges, 'FaceColor', clr1, 'EdgeColor', 'none', 'FaceAlpha', 1);
histogram(E0_filtered, bin_edges, 'FaceColor', clr2, 'EdgeColor', 'none', 'FaceAlpha', 1);
plot([temp_fit_unfiltered, temp_fit_unfiltered], [0, f1(f1.b)], 'Color', clr3, 'LineWidth', linw, 'LineStyle', '--')
plot([temp_fit_filtered, temp_fit_filtered], [0, f2(f2.b)], 'Color', clr4, 'LineWidth', linw, 'LineStyle', '--')
plot(bin_means, f1(bin_means), 'Color', clr3, 'LineWidth', linw)
plot(bin_means, f2(bin_means), 'Color', clr4, 'LineWidth', linw)
xlim([0, 3]); ylim([0, max([hist1, hist2])*1.05])
xlabel('Unaccelerated characteristic energy (keV)')
pbaspect(ar)
legend('Unfiltered', 'Filtered', ...
    sprintf('Unfiltered:%s%i eV', ' ', round(f1.b * 1e3, -1)), ...
    sprintf('Filtered:%s%i eV', ' ', round(f2.b * 1e3, -1)), ...
    'Location', 'northeast')
legend('boxoff')

set(gcf, 'Color', 'w', 'InvertHardcopy', 'off')
set([ax1, ax2, ax3], 'Color', clbg)

print(fig, 'plots/00_determining_srce', '-dpng', '-r96')
close all

im = imread('plots/00_determining_srce.png');
for i = 1:length(letter_pos)
    im = insertText(im, letter_pos(i, :), char(64 + i), ...
        'AnchorPoint', 'LeftTop', 'BoxOpacity', 0, 'FontSize', ftsz*2);
end
imwrite(im, 'plots/00_determining_srce.png', 'png')
