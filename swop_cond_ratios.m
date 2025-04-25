%%
scl.r = 1e+0; unt.r = '';       clm.r = 'D1A';
scl.s = 1e+0; unt.s = 'S';      clm.s = 'L18';
scl.n = 1e+0; unt.n = 'm^{-3}'; clm.n = 'L9';
scl.U = 1e3;  unt.U = 'mW/m^2'; clm.U = 'L19';
scl.x = 1e-3; unt.x = 'km';

lbl.x = sprintf('Mag. East (%s)', unt.x);
lbl.y = sprintf('Mag. North (%s)', unt.x);
lbl.z = sprintf('Mag. Up (%s)', unt.x);

%%
direc.A = fullfile('..', '..', 'public_html', 'Gemini3D', 'swop_20230210_35487_A_09_unacc');
direc.B = fullfile('..', '..', 'public_html', 'Gemini3D', 'swop_20230210_35487_A_09');
cfg.A = gemini3d.read.config(direc.A);
cfg.B = gemini3d.read.config(direc.B);
dat.A = gemini3d.read.frame(direc.A, 'time', cfg.A.times(end));
dat.B = gemini3d.read.frame(direc.B, 'time', cfg.B.times(end));
xg = gemini3d.read.grid(direc.A);

%%
x1 = xg.x1(3:end-2);
x2 = xg.x2(3:end-2);
x3 = xg.x3(3:end-2);
dx1 = xg.dx1h;
dx2 = xg.dx2h;
dx3 = xg.dx3h;
[X2, X3] = ndgrid(x2, x3);
[DX1, ~, ~] = ndgrid(dx1, dx2, dx3);

%%
for c = cell2mat(fieldnames(cfg)')
    sigP.(c) = h5read(dat.(c).filename, '/sigPall');
    sigH.(c) = -h5read(dat.(c).filename, '/sigHall');
    SIGP.(c) = squeeze(dot(sigP.(c), DX1));
    SIGH.(c) = squeeze(dot(sigH.(c), DX1));

    [E1, E2, E3] = gemscr.postprocess.pot2field(xg, dat.(c).Phitop);
    E0_fn = fullfile(direc.(c), cfg.(c).E0_dir, gemini3d.datelab(cfg.(c).times(end)) + ".h5");
    E2_BG = mean(h5read(E0_fn, '/Exit'), 'all');
    E3_BG = mean(h5read(E0_fn, '/Eyit'), 'all');
    E2 = E2 + E2_BG;
    E3 = E3 + E3_BG;
    j1 = dat.(c).J1;
    j2 = dat.(c).J2;
    j3 = dat.(c).J3;
    joule.(c) = (j1.*E1 + j2.*E2 + j3.*E3);
    joule_int.(c) = squeeze(dot(joule.(c), DX1));
end

%%
bound_A = h5read(fullfile(direc.A, 'ext', 'current.h5'), '/Boundary/Primary') * scl.x;
bound_B = h5read(fullfile(direc.A, 'ext', 'current.h5'), '/Boundary/Secondary') * scl.x;

%%
lim.x = [-1, 1] * 140;
lim.y = [-1, 1] * 85;
% lim.p = [quantile([SIGP.A(:); SIGP.B(:)], 0.01), ...
    % quantile([SIGP.A(:); SIGP.B(:)], 0.99)];
lim.p = [4.5, 20.5];
lim.h = [quantile([SIGH.A(:); SIGH.B(:)], 0.01), ...
    quantile([SIGH.A(:); SIGH.B(:)], 0.99)];
lim.U = [quantile([joule_int.A(:); joule_int.B(:)], 0.01), ...
    quantile([joule_int.A(:); joule_int.B(:)], 0.99)];
% lim.U = [0.8, 5.3] / scl.U;
lim.r = [0, 3];

close all
fntn = 'Arial';
fnts = 11 * 2;
clbg = [20, 21, 20] / 255;
cltx = [1, 1, 1];
reset(0)
set(0, 'defaultSurfaceEdgeColor', 'flat')
jules.tools.setall(0, 'FontName', fntn)
jules.tools.setall(0, 'FontSize', fnts)
jules.tools.setall(0, 'Multiplier', 1)
jules.tools.setall(0, 'XColor', cltx)
jules.tools.setall(0, 'YColor', cltx)
colorcet = @jules.tools.colorcet;

figure('Position', [70, 70, 1400, 700],'PaperPosition', [0, 0, 11.5, 7] * 2, ...
    'Color', clbg, 'InvertHardcopy', 'off')
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact')

% row 1
nexttile
hold on
pcolor(X2*scl.x, X3*scl.x, SIGP.A)
plot(bound_A(1, :), bound_A(2, :), 'k')
plot(bound_B(1, :), bound_B(2, :), '--k')
colormap(gca, colorcet(clm.s))
clb = colorbar;
clb.Location = 'northoutside';
clb.Label.String = sprintf('Pedersen Conductance (%s)', unt.s);
clb.Label.Color = cltx;
clb.Label.FontSize = fnts * 1.2;
clb.Label.FontWeight = 'bold';
xlim(lim.x); ylim(lim.y); clim(lim.p)
xlabel(lbl.x); ylabel([newline, newline, lbl.y])
text(-195, 0, 'Unaccelerated', 'Color', cltx, 'FontSize', fnts * 1.2, ...
    'FontWeight', 'bold', 'Rotation', 90, 'HorizontalAlignment', 'center')

nexttile
hold on
pcolor(X2*scl.x, X3*scl.x, SIGH.A ./ SIGP.A)
plot(bound_A(1, :), bound_A(2, :), 'k')
plot(bound_B(1, :), bound_B(2, :), '--k')
colormap(gca, colorcet(clm.r))
clb = colorbar;
clb.Location = 'northoutside';
clb.Label.String = sprintf('Hall / Pedersen Conductance Ratio');
clb.Label.Color = cltx;
clb.Label.FontSize = fnts * 1.2;
clb.Label.FontWeight = 'bold';
xlim(lim.x); ylim(lim.y); clim(lim.r)
xlabel(lbl.x); ylabel(lbl.y)

nexttile
hold on
pcolor(X2*scl.x, X3*scl.x, joule_int.A*scl.U)
plot(bound_A(1, :), bound_A(2, :), 'k')
plot(bound_B(1, :), bound_B(2, :), '--k')
colormap(gca, colorcet(clm.U))
clb = colorbar;
clb.Location = 'northoutside';
clb.Label.String = sprintf('Height-Integrated Joule Heating (%s)', unt.U);
clb.Label.Color = cltx;
clb.Label.FontSize = fnts * 1.2;
clb.Label.FontWeight = 'bold';
xlim(lim.x); ylim(lim.y); clim(lim.U*scl.U)
xlabel(lbl.x); ylabel(lbl.y)

% row 2
nexttile
hold on
pcolor(X2*scl.x, X3*scl.x, SIGP.B)
plot(bound_A(1, :), bound_A(2, :), 'k')
plot(bound_B(1, :), bound_B(2, :), '--k')
colormap(gca, colorcet(clm.s))
xlim(lim.x); ylim(lim.y); clim(lim.p)
xlabel(lbl.x); ylabel(lbl.y)
text(-195, 0, 'Accelerated (T_s = 490 eV)', 'Color', cltx, 'FontSize', fnts * 1.2, ...
    'FontWeight', 'bold', 'Rotation', 90, 'HorizontalAlignment', 'center')

nexttile
hold on
pcolor(X2*scl.x, X3*scl.x, SIGH.B ./ SIGP.B)
plot(bound_A(1, :), bound_A(2, :), 'k')
plot(bound_B(1, :), bound_B(2, :), '--k')
colormap(gca, colorcet(clm.r))
xlim(lim.x); ylim(lim.y); clim(lim.r)
xlabel(lbl.x); ylabel(lbl.y)

nexttile
hold on
pcolor(X2*scl.x, X3*scl.x, joule_int.B*scl.U)
plot(bound_A(1, :), bound_A(2, :), 'k')
plot(bound_B(1, :), bound_B(2, :), '--k')
colormap(gca, colorcet(clm.U))
xlim(lim.x); ylim(lim.y); clim(lim.U*scl.U)
xlabel(lbl.x); ylabel(lbl.y)

%%
filename = fullfile(direc.B, 'integrated_comparisons.png');
print(gcf, filename, '-dpng', '-r96')
close all

%%
p0 = [-30, -19, 111];
r0 = 4;
r1 = 4;
v0 = [0, 1, 0];
v1 = [0, 0, 1];
resolution = 200;

v0 = v0 / norm(v0);
v1 = v1 - dot(v0, v1) * v0;
v1 = v1 / norm(v1);
t = repmat(linspace(0, 1, resolution)', 1, 3);
c0 = p0 + r0*cos(2*pi*t).*v0 + r1*sin(2*pi*t).*v1;

[X1, XX3] = ndgrid(x1, x3);
[~, x2id] = min(abs(x2/1e3 - p0(1)));

close all
figure('Color',clbg)
tiledlayout(1, 2)

nexttile
hold on
pcolor(XX3*scl.x, X1*scl.x, squeeze(dat.A.J2(:, x2id, :))*1e6)
plot(c0(:, 2), c0(:, 3), 'r')
ylim([80, 250])
colorbar

nexttile
hold on
pcolor(XX3*scl.x, X1*scl.x, squeeze(dat.B.J2(:, x2id, :))*1e6)
plot(c0(:, 2), c0(:, 3), 'r')
ylim([80, 250])
colorbar