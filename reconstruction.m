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
lx2 = xg.lx(2); lx3 = xg.lx(3);
Bmag = double(abs(mean(xg.Bmag,'all')));

%%
time = datetime(ymd) + seconds(UTsec0 + recon_start);
dat = gemini3d.read.frame(direc,'time',time);
phi_tru = double(dat.Phitop);

%%
load(fullfile(direc,[orbit_fn,'_track.mat']),'tsat','x2sat','x3sat','v2sat','v3sat')
traj_ids = recon_start <= tsat & tsat <= recon_start + recon_dur;
sat_ids = 1:32;
% sat_ids = [1:4,29:32];

x2_traj = x2sat(traj_ids,sat_ids);
x3_traj = x3sat(traj_ids,sat_ids);
v2_traj = v2sat(traj_ids,sat_ids);
v3_traj = v3sat(traj_ids,sat_ids);

x2_traj = x2_traj(1:cad:end,:);
x3_traj = x3_traj(1:cad:end,:);
v2_traj = v2_traj(1:cad:end,:);
v3_traj = v3_traj(1:cad:end,:);

x_traj = [x2_traj(:),x3_traj(:)];
v_traj = [v2_traj(:),v3_traj(:)];

%%
prec_fn = gemini3d.datelab(time)+'.'+cfg.file_format;
prec_pth = fullfile(direc,cfg.prec_dir,prec_fn);
Q = double(h5read(prec_pth,'/Qp'));
Ebar = double(h5read(prec_pth,'/E0p'))/1e3;
SIGP = 40*Ebar.*sqrt(Q)./(16+Ebar.^2); % Robinson et al. (1987), Eq. (3)
edges = jules.tools.find_max_edges(SIGP);
x2_bound = x2(2:end-1);
x3_bound = nan(lx2-2,2);
for i = 1:lx2-2
    [~,bound_ids] = jules.tools.peak_detect(edges(i,:),num=2,smoothness=1);
    x3_bound(i,:) = x3(bound_ids);
end
bound = [x2_bound,x3_bound(:,2)];

%% reconstructing
phi_fit = jules.tools.reconstruct(x_traj,v_traj,bound,xg);

%% phi to E to v
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

%%
reg_bl = [x2_traj(1,1), x3_traj(1,1)]; % 4 corners of s/c array
reg_tl = [x2_traj(end,4), x3_traj(end,4)];
reg_br = [x2_traj(1,29), x3_traj(1,29)];
reg_tr = [x2_traj(end,32), x3_traj(end,32)];
regf_b = griddedInterpolant([reg_bl(1) reg_br(1)],[reg_bl(2) reg_br(2)]); % 4 boundaries of s/c array
regf_t = griddedInterpolant([reg_tl(1) reg_tr(1)],[reg_tl(2) reg_tr(2)]);
regf_l = griddedInterpolant([reg_bl(2) reg_tl(2)],[reg_bl(1) reg_tl(1)]);
regf_r = griddedInterpolant([reg_br(2) reg_tr(2)],[reg_br(1) reg_tr(1)]);
reg_b = false(lx2,lx3);
offset = 30e3;
for ix2 = 1:lx2
    for ix3 = 1:lx3
        x2p = X2(ix2,ix3);
        x3p = X3(ix2,ix3);
        if regf_l(x3p) - offset < x2p && x2p < regf_r(x3p) + offset % include coordinates left of right-most border, right of left-most, etc.
            if regf_b(x2p) - offset < x3p && x3p < regf_t(x2p) + offset
                reg_b(ix2,ix3) = true;
            end
        end
    end
end
reg_d = double(reg_b); % for plotting purposes
reg_d(reg_d==0) = nan;

v2_err_avg = mean(abs(v2_err(reg_b)));
v3_err_avg = mean(abs(v3_err(reg_b)));
v2_err_std = std(abs(v2_err(reg_b)));
v3_err_std = std(abs(v3_err(reg_b)));
v2_err_max = max(abs(v2_err(reg_b)));
v3_err_max = max(abs(v3_err(reg_b)));
fprintf('\n')
fprintf('Absolute eastward flow difference = %.0f +/- %.0f m/s\n',v2_err_avg,2*v2_err_std)
fprintf('Absolute northward flow difference = %.0f +/- %.0f m/s\n',v3_err_avg,2*v3_err_std)
fprintf('Maximum absolute eastward flow difference = %.0f m/s\n',v2_err_max)
fprintf('Maximum absolute northward flow difference = %.0f m/s\n',v3_err_max)

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
    jules.tools.hsv_params(V2_tru,V3_tru,MLAT,MLON,ALT,300e3,mlon_ref,hsv_sat);
[hsv_map_clb_fit,~,~,hsv_alt_fit,hsv_alt_map_fit] = ...
    jules.tools.hsv_params(V2_fit,V3_fit,MLAT,MLON,ALT,300e3,mlon_ref,hsv_sat);
[hsv_map_clb_err,~,~,hsv_alt_err,hsv_alt_map_err] = ...
    jules.tools.hsv_params(V2_err,V3_err,MLAT,MLON,ALT,300e3,mlon_ref,hsv_sat*0.3);

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
jules.tools.setall(0,'FontName',ftn)
jules.tools.setall(0,'FontSize',fts)
jules.tools.setall(0,'Multiplier',1)

colorcet = @jules.tools.colorcet;

scl.x = 1e-3; scl.v = 1e-3; scl.ve = 1e0; scl.p = 1e-3;
unt.x = 'km'; unt.v = 'km/s'; unt.ve = 'm/s'; unt.p = 'kV';
clm.v = 'D2'; clm.p = 'D10';

scl.vec = 50*scl.v;

lbl.x = sprintf('Mag. E (%s)',unt.x);
lbl.y = sprintf('Mag. N (%s)',unt.x);
lbl.p = sprintf('Potential (%s)',unt.p);
lbl.v = sprintf('v (%s)',unt.v);
lbl.ve = sprintf('\\Deltav (%s)',unt.ve);
lbl.vx = sprintf('v_E (%s)',unt.v);
lbl.vy = sprintf('v_N (%s)',unt.v);

lim.x = [min(x_traj(:,1)),max(x_traj(:,1))]*1.1*scl.x;
lim.y = [min(x_traj(:,2)),max(x_traj(:,2))]*1.4*scl.x;
lim.p = [min(phi_tru(:)),max(phi_tru(:))]*scl.p;
lim.vx = [-1,1]*max(abs(v_traj(:,1)))*scl.v;
lim.vy = [-1,1]*max(abs(v_traj(:,2)))*scl.v;
lim.ve = [-1,1]*max([v2_err_max,v3_err_max])*scl.ve;

ar = [range(lim.x),range(lim.y),range(lim.y)];

figure
set(gcf,'PaperPosition',[0,0,13.2,6])
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
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.vec,v3_traj*scl.vec,0,'.-r')
quiver(370,-205,1e3*scl.vec,0,0,'.-r')
text(430,-200,'1 km/s','FontSize',fts*0.8)
colormap(gca,hsv_alt_map_tru)
xlim(lim.x); ylim(lim.y); clim([0,1])
xticks([]); yticks([])
pbaspect(ar)

nexttile
title('Eastward flow difference')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,v2_err*scl.ve);
contour(X2*scl.x,X3*scl.x,reg_b,1,'Color',[1,1,1]*0.5)
plot(bound(:,1)*scl.x,bound(:,2)*scl.x,'--k')
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.ve)
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
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.vec,v3_traj*scl.vec,0,'.-r')
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
title('Northward flow difference')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,v3_err*scl.ve);
contour(X2*scl.x,X3*scl.x,reg_b,1,'Color',[1,1,1]*0.5)
plot(bound(:,1)*scl.x,bound(:,2)*scl.x,'--k')
colormap(gca,colorcet(clm.v))
clb = colorbar;
clb.Location = 'southoutside';
clb.Label.String = lbl.ve;
xlim(lim.x); ylim(lim.y); clim(lim.ve)
xlabel(lbl.x)
yticks([])
pbaspect(ar)

saveas(gcf,fullfile('plots','paper0','reconstruction.png'))
close all
