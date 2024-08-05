%% load simulation
direc = '\\dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\swop_0319_B_02';
cfg = gemini3d.read.config(direc);
time = cfg.times(end);
% xg = gemini3d.read.grid(direc);
dat = gemini3d.read.frame(direc,'time',time);

%% unpack grid
x2 = xg.x2(3:end-2);
x3 = xg.x3(3:end-2);
[X2,X3] = ndgrid(x2,x3);
lx1 = xg.lx(1); lx2 = xg.lx(2); lx3 = xg.lx(3);

mlon_to_x2 = griddedInterpolant(mlon,x2);
mlat_to_x3 = griddedInterpolant(mlat,x3);

% formatting simulation data
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
fSIGP_ext = griddedInterpolant(X2_ext,X3_ext,SIGP_ext,'linear');
fSIGH_ext = griddedInterpolant(X2_ext,X3_ext,SIGH_ext,'linear');
SIGP_ext = fSIGP_ext(X2,X3);
SIGH_ext = fSIGH_ext(X2,X3);

% add background electric fields
E0_fn = fullfile(direc,cfg.E0_dir,[char(gemini3d.datelab(time)),'.h5']);
E2_BG = mean(h5read(E0_fn,'/Exit'),'all');
E3_BG = mean(h5read(E0_fn,'/Eyit'),'all');
E2 = E2 + E2_BG;
E3 = E3 + E3_BG;

% implicit simulation data
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

%% plot
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

scl.x = 1e-3; scl.j = 1e6;
unt.j = 'uA/m^2';
clm.j = 'D1A';

jA_ext_p = jA_ext*scl.j; jB_ext_p = jB_ext*scl.j; jC_ext_p = jC_ext*scl.j;
jA_gem_p = jA_gem*scl.j; jB_gem_p = jB_gem*scl.j; jC_gem_p = jC_gem*scl.j;
deltaA = jA_gem_p - jA_ext_p;
deltaB = jB_gem_p - jB_ext_p;
deltaC = jC_gem_p - jC_ext_p;

lbl.x = 'Mag. E (km)';
lbl.y = 'Mag. N (km)';
lim.x = [-1,1]*max(x2)*scl.x*0.99;
lim.y = [-1,1]*max(x3)*scl.x*0.99;
lim.j = [-1,1]*quantile(abs([jA_gem_p(:);jB_gem_p(:);jC_gem_p(:)]),0.99);
ar = [1.618,1,1];

bound_h5 = fullfile(direc,'ext','potential.h5');
bound.A = h5read(bound_h5,'/Boundary/Primary');
bound.B = h5read(bound_h5,'/Boundary/Secondary');

figure
set(gcf,'PaperPosition',[0,0,10.2,14])
tlo = tiledlayout(3,1);
ltr = 'A';

% row 1
nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,jA_ext_p)
plot(bound.A(1,:)*scl.x,bound.A(2,:)*scl.x,'--k')
plot(bound.B(1,:)*scl.x,bound.B(2,:)*scl.x,'--k')
ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
xticks([])
clb = colorbar(gca);
clb.Label.String = sprintf('\\Sigma_P \\nabla_\\perp\\cdot\\bfE\\rm (%s)',unt.j);
colormap(gca,colorcet(clm.j))
clim(lim.j)
pbaspect(ar)

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,jB_ext_p)
plot(bound.A(1,:)*scl.x,bound.A(2,:)*scl.x,'--k')
plot(bound.B(1,:)*scl.x,bound.B(2,:)*scl.x,'--k')
ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
xticks([])
clb = colorbar(gca);
clb.Label.String = sprintf('\\bfE\\rm\\cdot\\nabla\\Sigma_P (%s)',unt.j);
colormap(gca,colorcet(clm.j))
clim(lim.j)
pbaspect(ar)

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,jC_ext_p)
plot(bound.A(1,:)*scl.x,bound.A(2,:)*scl.x,'--k')
plot(bound.B(1,:)*scl.x,bound.B(2,:)*scl.x,'--k')
xlabel(lbl.x); ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
clb = colorbar(gca);
clb.Label.String = sprintf('(\\bfE\\times b\\rm)\\cdot\\nabla\\Sigma_H (%s)',unt.j);
colormap(gca,colorcet(clm.j))
clim(lim.j)
pbaspect(ar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nexttile
% text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
% hold on
% pcolor(X2*scl.x,X3*scl.x,jA_gem_p)
% plot(bound.A(1,:)*scl.x,bound.A(2,:)*scl.x,'--k')
% plot(bound.B(1,:)*scl.x,bound.B(2,:)*scl.x,'--k')
% ylabel(lbl.y)
% xlim(lim.x); ylim(lim.y)
% xticks([])
% clb = colorbar(gca);
% clb.Label.String = sprintf('\\Sigma_P \\nabla_\\perp\\cdot\\bfE\\rm (%s)',unt.j);
% colormap(gca,colorcet(clm.j))
% clim(lim.j)
% pbaspect(ar)
% 
% nexttile
% text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
% hold on
% pcolor(X2*scl.x,X3*scl.x,jB_gem_p)
% plot(bound.A(1,:)*scl.x,bound.A(2,:)*scl.x,'--k')
% plot(bound.B(1,:)*scl.x,bound.B(2,:)*scl.x,'--k')
% ylabel(lbl.y)
% xlim(lim.x); ylim(lim.y)
% xticks([])
% clb = colorbar(gca);
% clb.Label.String = sprintf('\\bfE\\rm\\cdot\\nabla\\Sigma_P (%s)',unt.j);
% colormap(gca,colorcet(clm.j))
% clim(lim.j)
% pbaspect(ar)
% 
% nexttile
% text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
% hold on
% pcolor(X2*scl.x,X3*scl.x,jC_gem_p)
% plot(bound.A(1,:)*scl.x,bound.A(2,:)*scl.x,'--k')
% plot(bound.B(1,:)*scl.x,bound.B(2,:)*scl.x,'--k')
% xlabel(lbl.x); ylabel(lbl.y)
% xlim(lim.x); ylim(lim.y)
% clb = colorbar(gca);
% clb.Label.String = sprintf('(\\bfE\\times b\\rm)\\cdot\\nabla\\Sigma_H (%s)',unt.j);
% colormap(gca,colorcet(clm.j))
% clim(lim.j)
% pbaspect(ar)

[~,filename] = fileparts(direc);
saveas(gcf,fullfile('plots','swop',[filename,'_cond.png']))
close all
