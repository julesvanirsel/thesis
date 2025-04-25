glow_unacc = read_glow('data/dasc/23041_35487');
glow_acc = read_glow('data/dasc/23041_35487_acc');
ftsz = 15;
letter_pos = [[5, 0]; [800, 0]];

%%
Qp = 10; % mW/m2
E0 = 3; % keV
Ud = E0; % keV
Es = 0.490; % keV
maxalt = 190; % km

reset(0)
setall(0, 'FontName', 'Arial')
setall(0, 'FontSize', 10*2)
setall(0, 'Multiplier', 1)
setall(0, 'LineWidth', 1.5)
set(0, 'defaultAxesFontSizeMode', 'manual')
set(0, 'defaultSurfaceEdgeColor', 'flat')
set(0, 'defaultAxesLineWidth', 1)

close all
fig = figure;
tlo = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 6.5, 3]*2, ...
    'Position', [100, 100, 1200, 1200*3/6.5])

E = logspace(-1, 2, 1000);
phi_unacc = (1/(2*E0^2))*(E/E0).*exp(-E/E0);
phi_acc = (1/(Es^2+(Es+Ud)^2))*(E/Es).*exp(-(E-Ud)/Es);
phi_acc(E < Ud) = 0;

nexttile([1 2])
hold on
area(E, phi_unacc, 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r')
area(E, phi_acc, 'FaceColor', 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlim([0.2, 30]); ylim([1e-4, 5])
xlabel('Energy, E (keV)')
ylabel('Diff. Num. Flux, \phi/Q_p (keV^{-2})')
grid on
legend('Unaccelerated Maxwellian (E_0 = 3 keV)', ...
    'Accelerated Maxwellian (E_s = 490 eV, U_d = 3 keV)', ...
    'Location', 'northwest')

nexttile
hold on
clr = 'r';
for glow = [glow_unacc, glow_acc]
    [~, qid] = min(abs(glow.qvec - Qp));
    [~, eid] = min(abs(glow.ecvec/1e3 - E0));
    [~, altid] = min(abs(glow.altvec - maxalt));
    ne = squeeze(glow.ne(qid, eid, 1:altid));
    alt = glow.altvec(1:altid);

    plot(log10(ne), alt, clr)
    xlim([9.8, 11.8]); ylim([80, 189])
    clr = 'b';
end
xlabel('log_{10} Density, n_e (m^{-3})')
ylabel('Altitude (km)')
set(gca, 'YAxisLocation', 'right')
grid on

annotation('textarrow', [0.55, 0.51], [0.7, 0.54], ...
    'String', 'High energy tail', 'Color', 'r')
annotation('textarrow', [0.83, 0.89], [0.5, 0.39], ...
    'String', 'Lower density peak', 'Color', 'r')

print(fig, 'plots/00_spectra_comparison.png', '-dpng', '-r96');
close all

im = imread('plots/00_spectra_comparison.png');
for i = 1:length(letter_pos)
    im = insertText(im, letter_pos(i, :), char(64 + i), ...
        'AnchorPoint', 'LeftTop', 'BoxOpacity', 0, 'FontSize', ftsz*2);
end
imwrite(im, 'plots/00_spectra_comparison.png', 'png')