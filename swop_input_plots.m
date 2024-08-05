direc = '\\dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\swop_0210_A_02';
cfg = gemini3d.read.config(direc);
xg = gemini3d.read.grid(direc);
time = cfg.times(end);
dat = gemini3d.read.frame(direc,'time',time);

%%
mlon = squeeze(rad2deg(xg.phi(1,:,1)));
mlat = squeeze(90-rad2deg(xg.theta(1,1,:)));
x2 = xg.x2(3:end-2)/1e3;
x3 = xg.x3(3:end-2)/1e3;
[X2,X3] = ndgrid(x2,x3);

mlon_to_x2 = griddedInterpolant(mlon,x2);
mlat_to_x3 = griddedInterpolant(mlat,x3);

precip_h5 = fullfile(direc,'inputs','particles',[char(gemini3d.datelab(time)),'.h5']);
Qp = h5read(precip_h5,'/Qp');
Ep = h5read(precip_h5,'/E0p');

ext_h5 = fullfile(direc,'ext','precipitation.h5');
mlon_ext = h5read(ext_h5,'/Coordinates/Magnetic/Longitude');
mlat_ext = h5read(ext_h5,'/Coordinates/Magnetic/Latitude');
SIGP_ext = h5read(ext_h5,'/Derived/Conductance/Pedersen');
[X2_ext,X3_ext] = ndgrid(mlon_to_x2(mlon_ext),mlat_to_x3(mlat_ext));
fSIGP_ext = griddedInterpolant(X2_ext,X3_ext,SIGP_ext,'spline');
SIGP = fSIGP_ext(X2,X3);

bound_h5 = fullfile(direc,'ext','potential.h5');
bound.A = h5read(bound_h5,'/Boundary/Primary') / 1e3;
bound.B = h5read(bound_h5,'/Boundary/Secondary') / 1e3;

%%
close all
fts = 18*2; % fontsize
ftn = 'Arial';
lw = 2;

reset(0)
set(0,'defaultFigurePaperUnits','inches')
set(0,'defaultTiledlayoutPadding','tight')
set(0,'defaultTiledlayoutTileSpacing','compact')
set(0,'defaultLineLineWidth',lw)
jules.tools.setall(0,'FontName',ftn)
jules.tools.setall(0,'FontSize',fts)
jules.tools.setall(0,'Multiplier',1)
set(0,'defaultAxesFontSizeMode','manual')
set(0,'defaultSurfaceEdgeColor','flat')

colorcet = @jules.tools.colorcet;
scl.c = 1e-3; unt.c = 'keV';    clm.c = 'L17';
scl.U = 1e+3; unt.U = 'mW/m^2'; clm.U = 'L19';
scl.s = 1e+0; unt.s = 'S';      clm.s = 'L18';

lbl.x = 'Mag. E (km)';
lbl.y = 'Mag. N (km)';
lim.x = [-1,1]*max(x2)*0.99;
lim.y = [-1,1]*max(x3)*0.99;
ar = [1.618,1,1];

close all
figure('PaperPosition',[0,0,10,14])
tiledlayout(3,1)
ltr = 'A';

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2,X3,Qp)
plot(bound.A(1,:),bound.A(2,:),'--k')
plot(bound.B(1,:),bound.B(2,:),'--k')
clb = colorbar;
clb.Label.String = sprintf('Q_p (%s)',unt.U);
colormap(gca,colorcet(clm.U))
ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
xticks([])
pbaspect(ar)

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2,X3,Ep*scl.c)
plot(bound.A(1,:),bound.A(2,:),'--k')
plot(bound.B(1,:),bound.B(2,:),'--k')
clb = colorbar;
clb.Label.String = sprintf('E_p (%s)',unt.c);
colormap(gca,colorcet(clm.c))
ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
xticks([])
pbaspect(ar)

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2,X3,SIGP)
plot(bound.A(1,:),bound.A(2,:),'--k')
plot(bound.B(1,:),bound.B(2,:),'--k')
clb = colorbar;
clb.Label.String = sprintf('\\Sigma_P (%s)',unt.s);
colormap(gca,colorcet(clm.s))
xlabel(lbl.x); ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
pbaspect(ar)

[~,filename] = fileparts(direc);
saveas(gcf,fullfile('plots','swop',[filename,'_input.png']))
close all