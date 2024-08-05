%% load simulation
direc = '\\dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\swop_0210_AC_02';
cfg = gemini3d.read.config(direc);
time = cfg.times(end);
% xg = gemini3d.read.grid(direc);
dat = gemini3d.read.frame(direc,'time',time);

%% unpack grid
x2 = xg.x2(3:end-2);
x3 = xg.x3(3:end-2);
[X2,X3] = ndgrid(x2,x3);
lx1 = xg.lx(1); lx2 = xg.lx(2); lx3 = xg.lx(3);

% formatting simulation data
j1 = dat.J1;
phi = dat.Phitop;
[~,E2,E3] = gemscr.postprocess.pot2field(xg,phi);
[~,~,SIGP,SIGH] = jules.tools.load_conductances(direc,time,dat,cfg,xg);
SIGH = -SIGH;

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
[dSIGP3,dSIGP2] = gradient(SIGP);
[dSIGH3,dSIGH2] = gradient(SIGH);
jA = SIGP.*(dE22./dX2+dE33./dX3); % SIGP*div_perp(E)
jB = (dSIGP2.*squeeze(E2(1,:,:))./dX2 + dSIGP3.*squeeze(E3(1,:,:))./dX3); % grad(SIGP).E_perp
jC = (dSIGH2.*squeeze(E3(1,:,:))./dX2 - dSIGH3.*squeeze(E2(1,:,:))./dX3); % grad(SIGH).bxE_perp

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

jA_p = jA*scl.j; jB_p = jB*scl.j; jC_p = jC*scl.j;

lbl.x = 'Mag. E (km)';
lbl.y = 'Mag. N (km)';
lim.x = [-1,1]*max(x2)*scl.x*0.99;
lim.y = [-1,1]*max(x3)*scl.x*0.99;
lim.j = [-1,1]*quantile(abs([jA_p(:);jB_p(:);jC_p(:)]),0.99);
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
pcolor(X2*scl.x,X3*scl.x,jA_p)
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
pcolor(X2*scl.x,X3*scl.x,jB_p)
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
pcolor(X2*scl.x,X3*scl.x,0*jA_p+0*jB_p+jC_p)
plot(bound.A(1,:)*scl.x,bound.A(2,:)*scl.x,'--k')
plot(bound.B(1,:)*scl.x,bound.B(2,:)*scl.x,'--k')
xlabel(lbl.x); ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
clb = colorbar(gca);
clb.Label.String = sprintf('(\\bfE\\times b\\rm)\\cdot\\nabla\\Sigma_H (%s)',unt.j);
colormap(gca,colorcet(clm.j))
clim(lim.j)
pbaspect(ar)

[~,filename] = fileparts(direc);
saveas(gcf,fullfile('plots','swop',[filename,'_cc.png']))
close all
