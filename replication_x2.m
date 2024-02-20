if not(all(arrayfun(@exist,["in_situ","radar","image","xg","cfg"])))
    load('data\replicate_data_swop_03.mat')
end

%%
MLAT = 90-xg.theta*180/pi;
MLON = xg.phi*180/pi;
ALT = xg.alt;
x2 = double(xg.x2(3:end-2))';
x3 = double(xg.x3(3:end-2))';
[X2,X3] = ndgrid(x2,x3);
dx2 = xg.dx2h; dx3 = xg.dx3h;
lx2 = xg.lx(2); lx3 = xg.lx(3);
mlon_to_x2 = griddedInterpolant(squeeze(MLON(end,:,1)),x2);
mlat_to_x3 = griddedInterpolant(squeeze(MLAT(end,1,:)),x3);
x2_to_mlon = griddedInterpolant(x2,squeeze(MLON(end,:,1)));
x3_to_mlat = griddedInterpolant(x3,squeeze(MLAT(end,1,:)));
Bmag = abs(mean(xg.Bmag,'all'));
ar = [range(x2),range(x3),range(x3)];

clear('tracks')
tracks.swarm = in_situ;
tracks.pfisr = radar;
sig = 11;

%% replicate
[phi_both,~,~,E2_bg,E3_bg,~,~,weight,bound] = jules.tools.replicate(tracks,image,xg ...
    ,flow_smoothing_window = 10 ...
    ,boundary_smoothing_window = 30 ...
    ,show_plots = false ...
    ,save_plots = false ...
    ,direc = 'plots\paper0' ...
    ,suffix = 'track' ...
    ,add_phi_background = false ...
    ,fit_harmonic = true ...
    ,num_replications = 512 ...
    ,arc_definition = "conductance" ...
    ,edge_method = "contour" ...
    ,do_rotate = true ...
    ,do_scale = true ...
    ,contour_values = [1,1]*sig ...
    ,harmonic_mask = [1,1,1]*30e3 ...
    );

v2_bg = -E3_bg/Bmag;
v3_bg =  E2_bg/Bmag;

phi_swrm = jules.tools.replicate(tracks.swarm,image,xg ...
    ,flow_smoothing_window = 10 ...
    ,boundary_smoothing_window = 30 ...
    ,show_plots = false ...
    ,save_plots = false ...
    ,direc = 'plots\paper0' ...
    ,suffix = 'track' ...
    ,add_phi_background = false ...
    ,fit_harmonic = true ...
    ,num_replications = 512 ...
    ,arc_definition = "conductance" ...
    ,edge_method = "contour" ...
    ,do_rotate = true ...
    ,do_scale = true ...
    ,contour_values = [1,1]*sig ...
    ,harmonic_mask = [1,1,1]*30e3 ...
    ,flow_bg = [v2_bg,v3_bg] ...
    );

phi_pfsr = jules.tools.replicate(tracks.pfisr,image,xg ...
    ,flow_smoothing_window = 10 ...
    ,boundary_smoothing_window = 30 ...
    ,show_plots = false ...
    ,save_plots = false ...
    ,direc = 'plots\paper0' ...
    ,suffix = 'track' ...
    ,add_phi_background = false ...
    ,fit_harmonic = true ...
    ,num_replications = 512 ...
    ,arc_definition = "conductance" ...
    ,edge_method = "contour" ...
    ,do_rotate = true ...
    ,do_scale = true ...
    ,contour_values = [1,1]*sig ...
    ,harmonic_mask = [1,1,1]*30e3 ...
    ,flow_bg = [v2_bg,v3_bg] ...
    );

% phi to E to flow
if range(dx2) < 1e-3
    [E2_swrm,E3_swrm] = gradient(-phi_swrm',mean(dx2),mean(dx3));
    [E2_pfsr,E3_pfsr] = gradient(-phi_pfsr',mean(dx2),mean(dx3));
    [E2_both,E3_both] = gradient(-phi_both',mean(dx2),mean(dx3));
else
    [E2_swrm,E3_swrm] = gradient(-phi_swrm',dx2,dx3);
    [E2_pfsr,E3_pfsr] = gradient(-phi_pfsr',dx2,dx3);
    [E2_both,E3_both] = gradient(-phi_both',dx2,dx3);
end

v2_swrm = -E3_swrm'/Bmag;
v3_swrm =  E2_swrm'/Bmag;
v2_pfsr = -E3_pfsr'/Bmag;
v3_pfsr =  E2_pfsr'/Bmag;
v2_both = -E3_both'/Bmag;
v3_both =  E2_both'/Bmag;

%% plot
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
set(0,'defaultQuiverLineWidth',lw*0.1)
jules.tools.setall(0,'FontName',ftn)
jules.tools.setall(0,'FontSize',fts)
jules.tools.setall(0,'Multiplier',1)

colorcet = @jules.tools.colorcet;

scl.x = 1e-3; scl.v = 1e-3; scl.s = 1e0; scl.vec = 1e2; scl.q = 1e+3; 
unt.x = 'km'; unt.v = 'km/s'; unt.s = 'S'; unt.q = 'mW/m^2'; 
clm.v = 'D2'; clm.s = 'L18'; clm.q = 'L19';

lbl.x = sprintf('Mag. E (%s)',unt.x);
lbl.y = sprintf('Mag. N (%s)',unt.x);
lbl.vx = sprintf('v_E (%s)',unt.v);
lbl.vy = sprintf('v_N (%s)',unt.v);
lbl.s = sprintf('\\Sigma_P (%s)',unt.s);
lbl.q = sprintf('Q (%s)',unt.q);

qnt = 0.99;
max_v = quantile(abs([ ...
    [v2_swrm(:); v2_pfsr(:); v2_both(:)]+v2_bg ;...
    [v3_swrm(:); v3_pfsr(:); v3_both(:)]+v2_bg ...
    ]),qnt);
max_s = quantile(abs(image.pedersen(:)),qnt);
lim.v = [-1,1]*max_v*scl.v;
lim.s = [0,1]*max_s*scl.s;
lim.x = [-1,1]*max(x2)*scl.x;
lim.y = [-1,1]*max(x3)*scl.x;

% hsv plot variables
plot_bg = false;
hsv_sat = 700;
eps = 0.02;
mlon_ref = mean(MLON(:));

V2_swrm = permute(repmat(v2_swrm+v2_bg*plot_bg,[1,1,size(MLAT,1)]),[3,1,2]);
V3_swrm = permute(repmat(v3_swrm+v3_bg*plot_bg,[1,1,size(MLAT,1)]),[3,1,2]);
V2_pfsr = permute(repmat(v2_pfsr+v2_bg*plot_bg,[1,1,size(MLAT,1)]),[3,1,2]);
V3_pfsr = permute(repmat(v3_pfsr+v3_bg*plot_bg,[1,1,size(MLAT,1)]),[3,1,2]);
V2_both = permute(repmat(v2_both+v2_bg*plot_bg,[1,1,size(MLAT,1)]),[3,1,2]);
V3_both = permute(repmat(v3_both+v3_bg*plot_bg,[1,1,size(MLAT,1)]),[3,1,2]);
[hsv_map_clb,~,~,hsv_alt_swrm,hsv_alt_map_swrm] = ...
    jules.tools.hsv_params(V2_swrm,V3_swrm,MLAT,MLON,ALT,300e3,mlon_ref,hsv_sat);
[~,~,~,hsv_alt_pfsr,hsv_alt_map_pfsr] = ...
    jules.tools.hsv_params(V2_pfsr,V3_pfsr,MLAT,MLON,ALT,300e3,mlon_ref,hsv_sat);
[~,~,~,hsv_alt_both,hsv_alt_map_both] = ...
    jules.tools.hsv_params(V2_both,V3_both,MLAT,MLON,ALT,300e3,mlon_ref,hsv_sat);

X2_imag = mlon_to_x2(image.pos(:,:,1));
X3_imag = mlat_to_x3(image.pos(:,:,2));

%%
close all
figure
set(gcf,'PaperPosition',[0,0,13.2,5.6])
tiledlayout(2,3);
ltr = 'A';

x2_traj_swrm = mlon_to_x2(tracks.swarm.pos(:,1));
x3_traj_swrm = mlat_to_x3(tracks.swarm.pos(:,2));
v2_traj_swrm = tracks.swarm.flow(:,1) - v2_bg*(1-plot_bg);
v3_traj_swrm = tracks.swarm.flow(:,2) - v3_bg*(1-plot_bg);

x2_traj_pfsr = mlon_to_x2(tracks.pfisr.pos(:,1));
x3_traj_pfsr = mlat_to_x3(tracks.pfisr.pos(:,2));
v2_traj_pfsr = tracks.pfisr.flow(:,1) - v2_bg*(1-plot_bg);
v3_traj_pfsr = tracks.pfisr.flow(:,2) - v3_bg*(1-plot_bg);

x2_boundary = linspace(min(X2_imag(:)),max(X2_imag(:)),lx2);
x3_bounary_A = bound.A(x2_boundary);
x3_bounary_B = bound.B(x2_boundary);

% row 1
nexttile
title('Precip. flux')
text(0.04,1-0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2_imag*scl.x,X3_imag*scl.x,image.flux);
quiver(x2_traj_swrm*scl.x,x3_traj_swrm*scl.x, ...
    v2_traj_swrm*scl.v*scl.vec,v3_traj_swrm*scl.v*scl.vec,0,'.-b')
quiver(x2_traj_pfsr*scl.x,x3_traj_pfsr*scl.x, ...
    v2_traj_pfsr*scl.v*scl.vec,v3_traj_pfsr*scl.v*scl.vec,0,'.-b')
colormap(gca,colorcet(clm.q))
clb = colorbar;
clb.Label.String = lbl.q;
clb.Location = 'westoutside';
xlim(lim.x); ylim(lim.y);
ylabel(lbl.y)
xticks([])
pbaspect(ar)

nexttile
title('Pedersen cond.')
text(0.04,1-0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2_imag*scl.x,X3_imag*scl.x,image.pedersen);
quiver(x2_traj_swrm*scl.x,x3_traj_swrm*scl.x, ...
    v2_traj_swrm*scl.v*scl.vec,v3_traj_swrm*scl.v*scl.vec,0,'.-b')
quiver(x2_traj_pfsr*scl.x,x3_traj_pfsr*scl.x, ...
    v2_traj_pfsr*scl.v*scl.vec,v3_traj_pfsr*scl.v*scl.vec,0,'.-b')
plot(x2_boundary*scl.x,x3_bounary_A*scl.x,'k')
plot(x2_boundary*scl.x,x3_bounary_B*scl.x,'--k')
colormap(gca,colorcet(clm.s))
% clb = colorbar;
% clb.Label.String = lbl.s;
clim(lim.s)
xlim(lim.x); ylim(lim.y)
xticks([]); yticks([])
pbaspect(ar)

nexttile
title('Track A weighting')
text(0.04,1-0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,weight);
plot(x2_traj_swrm*scl.x,x3_traj_swrm*scl.x,'r','LineWidth',lw*0.1)
plot(x2_traj_pfsr*scl.x,x3_traj_pfsr*scl.x,'r','LineWidth',lw*0.1)
colormap(gca,colorcet('D3',reverse=true))
colorbar;
xlim(lim.x); ylim(lim.y)
xticks([]); yticks([])
pbaspect(ar)

% row 2
nexttile
title('Swarm flow')
text(0.04,1-0.9,char(ltr),'units','normalized','FontSize',fts*0.8,'Color','w'); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,hsv_alt_swrm);
quiver(x2_traj_swrm*scl.x,x3_traj_swrm*scl.x, ...
    v2_traj_swrm*scl.v*scl.vec,v3_traj_swrm*scl.v*scl.vec,0,'.-r')
plot(x2_traj_pfsr*scl.x,x3_traj_pfsr*scl.x,'r','LineWidth',0.1*lw)
contour(X2_imag*scl.x,X3_imag*scl.x,image.pedersen*scl.s,4,'--k')
colormap(gca,hsv_alt_map_swrm)
clim([0,1])
clb = colorbar;
colormap(clb,colorcet(clm.s))
clb.Label.String = lbl.s;
clb.Location = 'westoutside';
clb.Ticks = [0.1,0.5,0.9];
clb.TickLabels = round([0.1,0.5,0.9]*max_s);
xlim(lim.x); ylim(lim.y)
xlabel(lbl.x); ylabel(lbl.y)
pbaspect(ar)

nexttile
title('PFISR flow')
text(0.04,1-0.9,char(ltr),'units','normalized','FontSize',fts*0.8,'Color','w'); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,hsv_alt_pfsr);
plot(x2_traj_swrm*scl.x,x3_traj_swrm*scl.x,'r','LineWidth',lw*0.1)
quiver(x2_traj_pfsr*scl.x,x3_traj_pfsr*scl.x, ...
    v2_traj_pfsr*scl.v*scl.vec,v3_traj_pfsr*scl.v*scl.vec,0,'.-r')
contour(X2_imag*scl.x,X3_imag*scl.x,image.pedersen*scl.s,4,'--k')
colormap(gca,hsv_alt_map_pfsr)
clim([0,1])
xlim(lim.x); ylim(lim.y)
xlabel(lbl.x)
yticks([])
pbaspect(ar)

nexttile
title('Combined flow')
text(0.04,1-0.9,char(ltr),'units','normalized','FontSize',fts*0.8,'Color','w'); ltr = ltr +1;
hold on
pcolor(X2*scl.x,X3*scl.x,hsv_alt_both);
quiver(x2_traj_swrm*scl.x,x3_traj_swrm*scl.x, ...
    v2_traj_swrm*scl.v*scl.vec,v3_traj_swrm*scl.v*scl.vec,0,'.-r')
quiver(x2_traj_pfsr*scl.x,x3_traj_pfsr*scl.x, ...
    v2_traj_pfsr*scl.v*scl.vec,v3_traj_pfsr*scl.v*scl.vec,0,'.-r')
contour(X2_imag*scl.x,X3_imag*scl.x,image.pedersen*scl.s,4,'--k')
colormap(gca,hsv_alt_map_both)
clb = colorbar;
colormap(clb,hsv_map_clb)
clb.Limits = [0,1];
clb.Ticks = [0+eps,1/4,1/2,3/4,1-eps];
clb.TickLabels = {'W','S','E','N','W'};
clb.Label.String = ['Sat. at ',num2str(hsv_sat*scl.v),' ',unt.v];
clim([0,1])
xlim(lim.x); ylim(lim.y)
xlabel(lbl.x)
yticks([])
pbaspect(ar)

saveas(gcf,fullfile('plots','paper0','double_replication.png'))
close all