direc = '\\dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\swop_0210_AC_02';
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

%%
j1 = dat.J1;
phi = dat.Phitop;
[~,E2,E3] = gemscr.postprocess.pot2field(xg,phi);
[~,~,SIGP_gem,SIGH_gem] = jules.tools.load_conductances(direc,time,dat,cfg,xg);
SIGH_gem = -SIGH_gem;

%%
precip_h5 = fullfile(direc,'ext','precipitation.h5');

mlon_ext = h5read(precip_h5,'/Coordinates/Magnetic/Longitude');
mlat_ext = h5read(precip_h5,'/Coordinates/Magnetic/Latitude');
SIGP_ext = h5read(precip_h5,'/Derived/Conductance/Pedersen');
SIGH_ext = h5read(precip_h5,'/Derived/Conductance/Hall');

[X2_ext,X3_ext] = ndgrid(mlon_to_x2(mlon_ext),mlat_to_x3(mlat_ext));
fSIGP_ext = griddedInterpolant(X2_ext,X3_ext,SIGP_ext,'spline');
fSIGH_ext = griddedInterpolant(X2_ext,X3_ext,SIGH_ext,'spline');
SIGP_ext = fSIGP_ext(X2,X3);
SIGH_ext = fSIGH_ext(X2,X3);

%%
[~,dX2] = gradient(X2);
[dX3,~] = gradient(X3);
[~,dE22] = gradient(squeeze(E2(1,:,:)));
[dE33,~] = gradient(squeeze(E3(1,:,:)));

[dSIGP3_ext,dSIGP2_ext] = gradient(SIGP_ext);
[dSIGH3_ext,dSIGH2_ext] = gradient(SIGH_ext);
jA_ext = SIGP_ext.*(dE22./dX2+dE33./dX3); % SIGP*div_perp(E)
jB_ext = (dSIGP2_ext.*squeeze(E2(1,:,:))./dX2 + dSIGP3_ext.*squeeze(E3(1,:,:))./dX3); % grad(SIGP).E_perp
jC_ext = (dSIGH2_ext.*squeeze(E3(1,:,:))./dX2 - dSIGH3_ext.*squeeze(E2(1,:,:))./dX3); % grad(SIGH).bxE_perp

[dSIGP3_gem,dSIGP2_gem] = gradient(SIGP_gem);
[dSIGH3_gem,dSIGH2_gem] = gradient(SIGH_gem);
jA_gem = SIGP_gem.*(dE22./dX2+dE33./dX3); % SIGP*div_perp(E)
jB_gem = (dSIGP2_gem.*squeeze(E2(1,:,:))./dX2 + dSIGP3_gem.*squeeze(E3(1,:,:))./dX3); % grad(SIGP).E_perp
jC_gem = (dSIGH2_gem.*squeeze(E3(1,:,:))./dX2 - dSIGH3_gem.*squeeze(E2(1,:,:))./dX3); % grad(SIGH).bxE_perp

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
clm.j = 'D1A'; clm.s = 'L18';

lbl.x = 'Mag. E (km)';
lbl.y = 'Mag. N (km)';
lim.x = [-1,1]*max(x2)*0.99;
lim.y = [-1,1]*max(x3)*0.99;
% lim.j = [-1,1]*quantile(abs([jA_p(:);jB_p(:);jC_p(:)]),0.99);

close all
figure('PaperPosition',[0,0,10,14])
tiledlayout(3,1)
ltr = 'A';

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2,X3,dSIGP3_gem-dSIGP3_ext)
clb = colorbar;
clb.Label.String = sprintf('\\Delta(\\nabla_N\\Sigma_P) (S)');
colormap(gca,colorcet(clm.s))
ylabel(lbl.y)
xlim(lim.x); ylim(lim.y); clim([-1,1]*0.55)
xticks([])

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2,X3,dSIGH3_gem-dSIGH3_ext)
clb = colorbar;
clb.Label.String = sprintf('\\Delta(\\nabla_E\\Sigma_H) (S)');
colormap(gca,colorcet(clm.s))
ylabel(lbl.y)
xlim(lim.x); ylim(lim.y); clim([-1,1]*0.55)
xticks([])

[~,filename] = fileparts(direc);
saveas(gcf,fullfile('plots','swop',[filename,'_comp.png']))
close all
