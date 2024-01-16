%%
direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\aurora_longsharc_bent_01';
orbit_fn = 'arcs_orbit_20210211';
recon_start = 171;
recon_dur = 18;
cad = 32;

%%
cfg = gemini3d.read.config(direc);
ymd = cfg.ymd;
UTsec0 = cfg.UTsec0;
dtout = cfg.dtout;
tdur = cfg.tdur;

%%
xg = gemini3d.read.grid(direc);
x2 = double(xg.x2(3:end-2));
x3 = double(xg.x3(3:end-2));
dx2 = double(xg.dx2h);
dx3 = double(xg.dx3h);  
[X2,X3] = ndgrid(x2,x3);
lx2 = xg.lx(2);
Bmag = mean(xg.Bmag,'all');
ar = [range(x2),range(x3),range(x3)];

%%
time = datetime(ymd) + seconds(UTsec0 + recon_start);
dat = gemini3d.read.frame(direc,'time',time);
phi_tru = dat.Phitop;

%%
load(fullfile(direc,[orbit_fn,'_track.mat']),'tsat','x2sat','x3sat','v2sat','v3sat')
traj_ids = recon_start <= tsat & tsat <= recon_start + recon_dur;
x2_traj = x2sat(traj_ids,:);
x3_traj = x3sat(traj_ids,:);
v2_traj = v2sat(traj_ids,:);
v3_traj = v3sat(traj_ids,:);

x2_traj = x2_traj(1:cad:end,:);
x3_traj = x3_traj(1:cad:end,:);
v2_traj = v2_traj(1:cad:end,:);
v3_traj = v3_traj(1:cad:end,:);

x_traj = [x2_traj(:),x3_traj(:)];
v_traj = [v2_traj(:),v3_traj(:)];

%%
prec_fn = gemini3d.datelab(time)+'.'+cfg.file_format;
prec_pth = fullfile(direc,cfg.prec_dir,prec_fn);
Q = h5read(prec_pth,'/Qp');
Ebar = h5read(prec_pth,'/E0p')/1e3;
SIGP = 40*Ebar.*sqrt(Q)./(16+Ebar.^2); % Robinson et al. (1987), Eq. (3)
edges = tools.find_max_edges(SIGP);
x2_bound = x2(2:end-1);
x3_bound = nan(lx2-2,2);
for i = 1:lx2-2
    [~,bound_ids] = tools.peak_detect(edges(i,:),num=2,smoothness=1);
    x3_bound(i,:) = x3(bound_ids);
end
bound = [x2_bound,x3_bound(:,2)];

%%


%%
phi_fit = tools.reconstruct(x_traj,v_traj,bound,xg);

%%
pcolor(X2,X3,phi_fit); colorbar; shading flat

%%
[~,E2_tru,E3_tru] = gemscr.postprocess.pot2field(xg,phi_tru);
[~,E2_fit,E3_fit] = gemscr.postprocess.pot2field(xg,phi_fit);
v2_tru = -squeeze(E3_tru(end,:,:))/Bmag;
v3_tru = squeeze(E2_tru(end,:,:))/Bmag;
v2_fit = -squeeze(E3_fit(end,:,:))/Bmag;
v3_fit = squeeze(E2_fit(end,:,:))/Bmag;
v2_err = v2_fit-v2_tru;
v3_err = v3_fit-v3_tru;

fv2 = griddedInterpolant(X2,X3,v2_fit);
fv3 = griddedInterpolant(X2,X3,v3_fit);
v2_traj_fit = fv2(x2_traj,x3_traj);
v3_traj_fit = fv3(x2_traj,x3_traj);

%% hsv plot variables
MLAT = 90-squeeze(xg.theta)*180/pi;
MLON = squeeze(xg.phi)*180/pi;
ALT = xg.alt;
V2_tru = permute(repmat(v2_tru,[1,1,size(MLAT,1)]),[3,1,2]);
V3_tru = permute(repmat(v3_tru,[1,1,size(MLAT,1)]),[3,1,2]);
V2_fit = permute(repmat(v2_fit,[1,1,size(MLAT,1)]),[3,1,2]);
V3_fit = permute(repmat(v3_fit,[1,1,size(MLAT,1)]),[3,1,2]);
V2_err = permute(repmat(v2_err,[1,1,size(MLAT,1)]),[3,1,2]);
V3_err = permute(repmat(v3_err,[1,1,size(MLAT,1)]),[3,1,2]);
mlon_ref = mean(MLON(:));
hsv_sat = 700;
eps = 0.02;
[hsv_map_clb_tru,~,~,hsv_alt_tru,hsv_alt_map_tru] = ...
    tools.hsv_params(V2_fit,V3_fit,MLAT,MLON,ALT,300e3,mlon_ref,hsv_sat);
[hsv_map_clb_fit,~,~,hsv_alt_fit,hsv_alt_map_fit] = ...
    tools.hsv_params(V2_fit,V3_fit,MLAT,MLON,ALT,300e3,mlon_ref,hsv_sat);
[hsv_map_clb_err,~,~,hsv_alt_err,hsv_alt_map_err] = ...
    tools.hsv_params(V2_err,V3_err,MLAT,MLON,ALT,300e3,mlon_ref,hsv_sat*0.3);

%%
ftn = 'Arial';
fts = 10*2;
lw = 1.4;

close all
reset(0)
set(0,'defaultFigurePaperUnits','inches')
set(0,'defaultTiledlayoutPadding','tight')
set(0,'defaultTiledlayoutTileSpacing','tight')
set(0,'defaultSurfaceEdgeColor','flat')
set(0,'defaultLineLineWidth',lw)
set(0,'defaultScatterLineWidth',lw)
set(0,'defaultQuiverLineWidth',lw*0.7)
tools.setall(0,'FontName',ftn)
tools.setall(0,'FontSize',fts)
tools.setall(0,'Multiplier',1)

scl.x = 1e-3; scl.v = 1e-3; scl.p = 1e-3;
unt.x = 'km'; unt.v = 'km/s'; unt.p = 'kV';
clm.v = 'D2'; clm.p = 'D10';

lbl.x = sprintf('M. east (%s)',unt.x);
lbl.y = sprintf('M. north (%s)',unt.x);
lbl.p = sprintf('Potential (%s)',unt.p);
lbl.v = sprintf('v (%s)',unt.v);
lbl.vx = sprintf('v_{east} (%s)',unt.v);
lbl.vy = sprintf('v_{north} (%s)',unt.v);

lim.x = [min(x_traj(:,1)),max(x_traj(:,1))]*1.3*scl.x;
lim.y = [min(x_traj(:,2)),max(x_traj(:,2))]*1.3*scl.x;
lim.p = [min(phi_tru(:)),max(phi_tru(:))]*scl.p;
lim.vx = [-1,1]*max(abs(v_traj(:,1)))*scl.v;
lim.vy = [-1,1]*max(abs(v_traj(:,2)))*scl.v;

sc = 50;

figure
set(gcf,'PaperPosition',[0,0,13.2,5.5])
tiledlayout(2,3)
ltr = 65;

% row 1
nexttile
title('OSSE potential')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,phi_tru*scl.p)
plot(bound(:,1)*scl.x,bound(:,2)*scl.x,'--k')
colormap(gca,colorcet(clm.p))
xlim(lim.x); ylim(lim.y); clim(lim.p)
ylabel(lbl.y)
xticks([])
pbaspect(ar)

nexttile
title('OSSE flow')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,hsv_alt_tru);
plot(bound(:,1)*scl.x,bound(:,2)*scl.x,'--k')
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,hsv_alt_map_tru)
xlim(lim.x); ylim(lim.y); clim([0,1])
xticks([]); yticks([])
pbaspect(ar)

nexttile
title('Eastward flow error')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,v2_err*scl.v);
plot(bound(:,1)*scl.x,bound(:,2)*scl.x,'--k')
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.vx/5)
xticks([]); yticks([])
pbaspect(ar)

% row 2
nexttile
title('Reconstructed potential')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,phi_fit*scl.p)
plot(bound(:,1)*scl.x,bound(:,2)*scl.x,'--k')
colormap(gca,colorcet(clm.p))
clb = colorbar;
clb.Location = 'southoutside';
clb.Label.String = lbl.p;
xlim(lim.x); ylim(lim.y); clim(lim.p)
xlabel(lbl.x); ylabel(lbl.y)
pbaspect(ar)

nexttile
title('Reconstructed flow')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,hsv_alt_fit);
plot(bound(:,1)*scl.x,bound(:,2)*scl.x,'--k')
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,hsv_alt_map_fit)
clb = colorbar;
colormap(clb,hsv_map_clb_fit)
clb.Limits = [0,1];
clb.Ticks = [0+eps,1/4,1/2,3/4,1-eps];
clb.TickLabels = {'W','S','E','N','W'};
clb.Label.String = ['Sat. at ',num2str(hsv_sat*scl.v),' ',unt.v];
clb.Location = 'southoutside';
xlim(lim.x); ylim(lim.y); clim([0,1])
xlabel(lbl.x)
yticks([])
pbaspect(ar)

nexttile
title('Northward flow error')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,v3_err*scl.v);
plot(bound(:,1)*scl.x,bound(:,2)*scl.x,'--k')
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
clb = colorbar;
clb.Location = 'southoutside';
clb.Label.String = sprintf('\\Delta %s',lbl.v);
xlim(lim.x); ylim(lim.y); clim(lim.vx/5)
xlabel(lbl.x)
yticks([])
pbaspect(ar)

saveas(gcf,fullfile('plots','reconstruction.png'))
close all
