%% load simulation data - mar 19, 2023 - 8:23:44.2 UT
close all
direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\swop_03';
cfg = gemini3d.read.config(direc);

% unpack grid
xg = gemini3d.grid.cartesian(cfg);
mlon = squeeze(xg.phi(1,:,1)*180/pi);
mlat = squeeze(90-xg.theta(1,1,:)*180/pi);
x2 = double(xg.x2(3:end-2));
x3 = double(xg.x3(3:end-2));
dx2 = double(xg.dx2h);
dx3 = double(xg.dx3h);
mlon_to_x2 = griddedInterpolant(mlon,x2);
mlat_to_x3 = griddedInterpolant(mlat,x3);
Bmag = mean(xg.Bmag(:));

colorcet = @jules.tools.colorcet;

%% pack up image
load('data\SW_EXPT_EFIB_TCT02_20230319T081600_sky_inversion.mat','event')

lims_mlon = [min(mlon),max(mlon)];
lims_mlat = [min(mlat),max(mlat)];
lims_mlon = [min(event.maglon_110(:))+360,max(event.maglon_110(:))+360];
lims_mlat = [min(event.maglat_110(:)),max(event.maglat_110(:))];
% lims_mlon = [263,276];
% lims_mlat = [65.5,69];
% lims_mlon = [0,360];
% lims_mlat = [0,90];
[~,ida1] = min(abs(event.maglon_110(1,:)+360-lims_mlon(1)));
[~,ida2] = min(abs(event.maglon_110(1,:)+360-lims_mlon(2)));
[~,idb1] = min(abs(event.maglat_110(:,1)-lims_mlat(1)));
[~,idb2] = min(abs(event.maglat_110(:,1)-lims_mlat(2)));

clear('image')
image.pos(:,:,1) = event.maglon_110(idb1:idb2,ida1:ida2)'+360;
image.pos(:,:,2) = event.maglat_110(idb1:idb2,ida1:ida2)';
image.flux = event.Q(idb1:idb2,ida1:ida2)';
image.energy = event.E0(idb1:idb2,ida1:ida2)';
image.pedersen = event.SigP(idb1:idb2,ida1:ida2)';

% pack up in_situ
clear('in_situ')
load('data\SW_EXPT_EFIB_TCT02_20230319T081600_alongtrackincluded_UNTRUSTWORTHY.mat' ...
    ,'SW_EXPT_EFIB_TCT02_20230319T081600')
SW_EXPT_B = SW_EXPT_EFIB_TCT02_20230319T081600;

mlat_buf = 1;
ok = SW_EXPT_B.Mag_Lat > (lims_mlat(1)-mlat_buf) & ...
    SW_EXPT_B.Mag_Lat < (lims_mlat(2)+mlat_buf);
in_situ.pos(:,1) = SW_EXPT_B.Mag_Lon(ok) + 360;
in_situ.pos(:,2) = SW_EXPT_B.Mag_Lat(ok);
in_situ.flow(:,1) = SW_EXPT_B.Flow_MagEast(ok);
in_situ.flow(:,2) = SW_EXPT_B.Flow_MagNorth(ok);

% pack up radar
clear('radar')
path = 'data\20230319.001_lp_3min-fitcal-vvels_lat.h5';
times = h5read(path,'/Time/UnixTime');
tid = 34;
% tid = 21;
fprintf('Start time = %s\n',datetime(times(1,tid),'ConvertFrom','posixtime'))
fprintf('End time = %s\n',datetime(times(2,tid),'ConvertFrom','posixtime'))
mlon_vvels = h5read(path,'/VvelsMagCoords/Longitude')+360;
mlat_vvels = h5read(path,'/VvelsMagCoords/Latitude');
flow_vvels = h5read(path,'/VvelsMagCoords/Velocity');
flow_vvels_err = h5read(path,'/VvelsMagCoords/errVelocity');

ok1 = mean(flow_vvels_err(:,:,tid))' < 300;
ok2 = mlat_vvels > (lims_mlat(1)-mlat_buf) ...
    & mlat_vvels < (lims_mlat(2)+mlat_buf);
ok = ok1 & ok2;
l_vvels = length(mlon_vvels(ok));
radar.pos(:,1) = mlon_vvels(ok) + 1e-5*(1:l_vvels)';
radar.pos(:,2) = mlat_vvels(ok);
radar.flow(:,1) = squeeze(flow_vvels(1,ok,tid))';
radar.flow(:,2) = -squeeze(flow_vvels(2,ok,tid))';

% save('data\replicate_data_swop_03.mat','image','in_situ','radar','xg','cfg')

% check region
disp(any(isnan(image.flux(:))))
disp(any(isnan(image.energy(:))))
disp(any(isnan(image.pedersen(:))))

% plot
[~,sid] = min(abs(in_situ.pos(:,2)-66.8337));
sc = 0.002;
sig = 11;
ar = [range(x2),range(x3),1];
M = contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,[1,1]*sig);
[xa,ya] = jules.tools.get_contour(M,rank=1);
[xb,yb] = jules.tools.get_contour(M,rank=2);
close all

figure
set(gcf,'PaperPosition',[0,0,6.5,5.5])
set(gcf,'Position',[30,60,1500,800]);
set(gca,'Fontsize',20)
tiledlayout(2,2)

nexttile
hold on
pcolor(image.pos(:,:,1),image.pos(:,:,2),image.flux)
quiver(in_situ.pos(:,1),in_situ.pos(:,2),in_situ.flow(:,1)*sc,in_situ.flow(:,2)*sc,0,'.-k')
quiver(radar.pos(:,1),radar.pos(:,2),radar.flow(:,1)*sc,radar.flow(:,2)*sc,0,'.-b')
scatter(in_situ.pos(sid,1),in_situ.pos(sid,2),100,'xm')
xlim(lims_mlon); ylim(lims_mlat)
xlabel('M. Longitude'); ylabel('M. Latitude')
clb = colorbar;
clb.Label.String = 'Q (mW/m^2)';
clim([0,32])
colormap(gca,colorcet('L18'))
pbaspect(ar)
shading flat

nexttile
hold on
pcolor(image.pos(:,:,1),image.pos(:,:,2),image.energy/1e3)
quiver(in_situ.pos(:,1),in_situ.pos(:,2),in_situ.flow(:,1)*sc,in_situ.flow(:,2)*sc,0,'.-k')
quiver(radar.pos(:,1),radar.pos(:,2),radar.flow(:,1)*sc,radar.flow(:,2)*sc,0,'.-b')
scatter(in_situ.pos(sid,1),in_situ.pos(sid,2),100,'xm')
xlim(lims_mlon); ylim(lims_mlat)
xlabel('M. Longitude'); ylabel('M. Latitude')
clb = colorbar;
clb.Label.String = 'E_0 (keV)';
clim([0,6.5])
colormap(gca,colorcet('L17'))
pbaspect(ar)
shading flat

nexttile
hold on
pcolor(image.pos(:,:,1),image.pos(:,:,2),image.flux);
contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,10)
quiver(in_situ.pos(:,1),in_situ.pos(:,2),in_situ.flow(:,1)*sc,in_situ.flow(:,2)*sc,0,'.-k')
quiver(radar.pos(:,1),radar.pos(:,2),radar.flow(:,1)*sc,radar.flow(:,2)*sc,0,'.-b')
scatter(in_situ.pos(sid,1),in_situ.pos(sid,2),100,'xm')
xlim(lims_mlon); ylim(lims_mlat)
xlabel('M. Longitude'); ylabel('M. Latitude')
clb = colorbar;
clb.Label.String = 'Q (mW/m^2)  and  \Sigma_P (S)';
clim([0,32])
colormap(gca,colorcet('L18'))
pbaspect(ar)
shading flat

nexttile
hold on
pcolor(image.pos(:,:,1),image.pos(:,:,2),image.flux);
plot(xa,ya,'k')
plot(xb,yb,'g')
contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,[1,1]*sig,'--m');
quiver(in_situ.pos(:,1),in_situ.pos(:,2),in_situ.flow(:,1)*sc,in_situ.flow(:,2)*sc,0,'.-k')
quiver(radar.pos(:,1),radar.pos(:,2),radar.flow(:,1)*sc,radar.flow(:,2)*sc,0,'.-b')
scatter(in_situ.pos(sid,1),in_situ.pos(sid,2),100,'xm')
xlim(lims_mlon); ylim(lims_mlat)
xlabel('M. Longitude'); ylabel('M. Latitude')
clb = colorbar;
clb.Label.String = 'Q (mW/m^2)  and  \Sigma_P (S)';
clim([0,32])
colormap(gca,colorcet('L18'))
pbaspect(ar)
shading flat

%% replicate
[phi_rep,mlon,mlat,E2_bg,E3_bg] = jules.tools.replicate(in_situ,image,xg ...
    ,flow_smoothing_window = 1 ...
    ,boundary_smoothing_window = 30 ...
    ,show_plots = false ...
    ,save_plots = false ...
    ,direc = 'plots' ...
    ,add_phi_background = false ...
    ,fit_harmonic = true ...
    ,num_replications = 512 ...
    ,arc_definition = "conductance" ...
    ,edge_method = "contour" ...
    ,do_rotate = true ...
    ,do_scale = true ...
    ,contour_values = [1,1]*sig ...
    ,harmonic_mask = [1,0,1]*50e3 ...
    ...%,flow_bg = [300,-100] ...
    );

%% reconstruct
% pos = in_situ.pos;
% flow = in_situ.flow;
% pos = [in_situ.pos; radar.pos];
% flow = [in_situ.flow; radar.flow];

load('data\reps.mat','x2_traj_rep','x3_traj_rep','v2_traj_rep','v3_traj_rep')
cd = 64;
x2_traj_rep = x2_traj_rep(cd/2:cd:end,:);
x3_traj_rep = x3_traj_rep(cd/2:cd:end,:);
v2_traj_rep = v2_traj_rep(cd/2:cd:end,:);
v3_traj_rep = v3_traj_rep(cd/2:cd:end,:);
pos = [x2_traj_rep(:), x3_traj_rep(:)];
flow = [v2_traj_rep(:), v3_traj_rep(:)];
quiver(pos(:,1),pos(:,2),flow(:,1),flow(:,2))

figure(99)
M = contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,[1,1]*sig);
close(gcf)
[bx,by] = jules.tools.get_contour(M,rank=2);
boundary = [bx; by]';

%%
[phi_rec,P] = jules.tools.reconstruct(pos,flow,boundary,xg);

%% plot
if range(dx2) < 1e-3
    [E2_rep,E3_rep] = gradient(-phi_rep',mean(dx2),mean(dx3));
    [E2_rec,E3_rec] = gradient(-phi_rec',mean(dx2),mean(dx3));
else
    [E2_rep,E3_rep] = gradient(-phi_rep',dx2,dx3);
    [E2_rec,E3_rec] = gradient(-phi_rec',dx2,dx3);
end
E2_rec = E2_rec - 0*E2_bg;
E3_rec = E3_rec - 0*E3_bg;

v2_rep = -E3_rep'/Bmag;
v3_rep =  E2_rep'/Bmag;
v2_rec = -E3_rec'/Bmag;
v3_rec =  E2_rec'/Bmag;
 
qnt = 0.99;
v2_max = quantile(abs([v2_rep(:);v2_rec(:)]),qnt);
v3_max = quantile(abs([v3_rep(:);v3_rec(:)]),qnt);
scl.x = 1e-3; scl.v = 1e-3;
lim.v2 = [-1,1]*v2_max*scl.v;
lim.v3 = [-1,1]*v3_max*scl.v;
close all

figure
set(gcf,'PaperPosition',[0,0,6.5,5.5])
set(gcf,'Position',[0,60,1920,1080-150]);
set(gca,'Fontsize',20)
tiledlayout(2,2)

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,v2_rep*scl.v)
clim(lim.v2)
pbaspect(ar)

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,v2_rec*scl.v)
colorbar;
clim(lim.v2)
pbaspect(ar)

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,v3_rep*scl.v)
clim(lim.v2)
pbaspect(ar)

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,v3_rec*scl.v)
colorbar;
clim(lim.v2)
pbaspect(ar)


%% old
% load('data\replicate_data_isinglass.mat','image','in_situ')
% image_IG = image;
% in_situ_IG = in_situ;
% clear('image','in_situ')

% load('data\0212_precip_map_lowres.mat')
% load('data\0314_precip_map.mat')
% load('data\SW_EXPT_EFIB_TCT02_20230319T081600_sky_inversion.mat','event')
% event = dt_20230314T0649;

% lims_mlon = [268.84, 276.26];
% lims_mlat = [63.65, 67.70];
% lims_mlon = [min(mlon),max(mlon)];
% lims_mlat = [min(mlat),max(mlat)];
% lims_mlon = [min(event.maglon_110(:))+360,max(event.maglon_110(:))+360];
% lims_mlat = [min(event.maglat_110(:)),max(event.maglat_110(:))];

% [~,ida1] = min(abs(event.maglon_110(1,:)+360-lims_mlon(1)));
% [~,ida2] = min(abs(event.maglon_110(1,:)+360-lims_mlon(2)));
% [~,idb1] = min(abs(event.maglat_110(:,1)-lims_mlat(1)));
% [~,idb2] = min(abs(event.maglat_110(:,1)-lims_mlat(2)));
% 
% clear('image')
% image.pos(:,:,1) = event.maglon_110(idb1:idb2,ida1:ida2)'+360;
% image.pos(:,:,2) = event.maglat_110(idb1:idb2,ida1:ida2)';
% image.flux = event.Q(idb1:idb2,ida1:ida2)';
% image.energy = event.E0(idb1:idb2,ida1:ida2)';
% image.pedersen = event.SigP(idb1:idb2,ida1:ida2)';
% clear('event')

% make in_situ
% load('data\SW_EXPT_EFIA_TCT02_20230314T064200.mat')
% load('data\SW_EXPT_EFIC_TCT02_20230314T064200.mat')
% load('data\SW_EXPT_EFIC_TCT02_20230212T101600.mat')
% load('data\SW_EXPT_EFIB_TCT02_20230319T081600_alongtrackincluded_UNTRUSTWORTHY.mat')
% load('data\SW_EXPT_EFIB_TCT02_20230319T081600_crosstrackonly.mat')
% SW_EXPT_B = SW_EXPT_EFIB_TCT02_20230319T081600;

% SW_EXPT_A = SW_EXPT_EFIA_TCT02_20230314T064200;
% SW_EXPT_C = SW_EXPT_EFIC_TCT02_20230314T064200;
% SW_EXPT_C = SW_EXPT_EFIC_TCT02_20230212T101600;

% in_situ_A.pos(:,1) = SW_EXPT_A.Mag_Lon + 360;
% in_situ_A.pos(:,2) = SW_EXPT_A.Mag_Lat;
% in_situ_A.flow(:,1) = SW_EXPT_A.Flow_MagEast;
% in_situ_A.flow(:,2) = SW_EXPT_A.Flow_MagNorth;
% clear('SW_EXPT_A')

% in_situ_B.pos(:,1) = SW_EXPT_B.Mag_Lon + 360;
% in_situ_B.pos(:,2) = SW_EXPT_B.Mag_Lat;
% in_situ_B.flow(:,1) = SW_EXPT_B.Flow_MagEast;
% in_situ_B.flow(:,2) = SW_EXPT_B.Flow_MagNorth;
% clear('SW_EXPT_B')

% in_situ_C.pos(:,1) = SW_EXPT_C.Mag_Lon+360;
% in_situ_C.pos(:,2) = SW_EXPT_C.Mag_Lat;
% in_situ_C.flow(:,1) = SW_EXPT_C.Flow_MagEast;
% in_situ_C.flow(:,2) = SW_EXPT_C.Flow_MagNorth;
% clear('SW_EXPT_C')

% save('data\replicate_data_swop_01.mat','image','in_situ_C')
% save('data\replicate_data_swop_03.mat','image','in_situ_B')

% fprintf('Altitude = %i m\n',galt(1,aid));
% 
% glon = squeeze(glon(:,aid))';
% glat = squeeze(glat(:,aid))';
% gv_east = squeeze(gv(1,:,aid,tid))/1e3; % see Laundal & Richmond (2016), Eq. 57-59
% gv_north = squeeze(gv(2,:,aid,tid))/1e3;
% gv_up = squeeze(gv(3,:,aid,tid))/1e3;

%
% close all
% 
% figure
% sid = 928;
% hold on
% pcolor(image.pos(:,:,1),image.pos(:,:,2),image.flux); colorbar
% % quiver(in_situ_A.pos(:,1),in_situ_A.pos(:,2),in_situ_A.flow(:,1),in_situ_A.flow(:,2),'.-r')
% % quiver(in_situ_C.pos(:,1),in_situ_C.pos(:,2),in_situ_C.flow(:,1),in_situ_C.flow(:,2),'.-r')
% quiver(in_situ_B.pos(:,1),in_situ_B.pos(:,2),in_situ_B.flow(:,1)*sc,in_situ_B.flow(:,2)*sc,0,'.-r')
% quiver(mlon_vvels,mlat_vvels,flow_east*sc,flow_north*sc,0,'.-y')
% scatter(in_situ_B.pos(sid,1),in_situ_B.pos(sid,2),100,'xm')
% xlim(lims_mlon); ylim(lims_mlat)
% pbaspect([range(x2),range(x3),1])
% shading flat
% 
% figure
% hold on
% pcolor(image.pos(:,:,1),image.pos(:,:,2),image.energy); colorbar
% % quiver(in_situ_A.pos(:,1),in_situ_A.pos(:,2),in_situ_A.flow(:,1),in_situ_A.flow(:,2),'.-r')
% % quiver(in_situ_C.pos(:,1),in_situ_C.pos(:,2),in_situ_C.flow(:,1),in_situ_C.flow(:,2),'.-r')
% quiver(in_situ_B.pos(:,1),in_situ_B.pos(:,2),in_situ_B.flow(:,1)*sc,in_situ_B.flow(:,2)*sc,0,'.-r')
% quiver(mlon_vvels,mlat_vvels,flow_east*sc,flow_north*sc,0,'.-y')
% scatter(in_situ_B.pos(sid,1),in_situ_B.pos(sid,2),100,'xm')
% xlim(lims_mlon); ylim(lims_mlat)
% pbaspect([range(x2),range(x3),1])
% shading flat

%
% figure
% set(gcf,'PaperPosition',[0,0,6.5,5.5])
% set(gca,'Fontsize',20)
% sc = 0.001;
% sc2 = 0.0005;
% % ids = 1:3:length(in_situ_A.flow);
% ids = 1:3:length(in_situ_B.flow);
% 
% hold on
% pcolor(image.pos(:,:,1),image.pos(:,:,2),image.flux); shading flat; colormap(colorcet('L18'))
% contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,10)
% % quiver(in_situ_A.pos(ids,1),in_situ_A.pos(ids,2),in_situ_A.flow(ids,1)*sc,in_situ_A.flow(ids,2)*sc,0,'.-b')
% % quiver(in_situ_C.pos(ids,1),in_situ_C.pos(ids,2),in_situ_C.flow(ids,1)*sc2,in_situ_C.flow(ids,2)*sc2,0,'.-b')
% quiver(in_situ_B.pos(ids,1),in_situ_B.pos(ids,2),in_situ_B.flow(ids,1)*sc,in_situ_B.flow(ids,2)*sc,0,'.-b')
% quiver(mlon_vvels,mlat_vvels,flow_east*sc,flow_north*sc,0,'.-k')
% xlim(lims_mlon); ylim(lims_mlat)
% xlabel('M. Longitude'); ylabel('M. Latitude')
% clb = colorbar;
% clb.Label.String = 'Q (mW/m^2)  and  \Sigma_P (S)';
% clim([0,30])
% pbaspect([1,1,1])
% 
% saveas(gcf,'plots/swop_dasc.png')

%
% disp(any(isnan(image.flux(:))))
% disp(any(isnan(image.energy(:))))
% disp(any(isnan(image.pedersen(:))))

%
% figure
% hold on
% contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,20); colorbar; shading flat
% contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,[1,1]*13,'m'); colorbar; shading flat
% % contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,[1,1]*13,'k'); colorbar; shading flat
% % quiver(in_situ_A.pos(:,1),in_situ_A.pos(:,2),in_situ_A.flow(:,1),in_situ_A.flow(:,2),'.-r')
% % quiver(in_situ_C.pos(:,1),in_situ_C.pos(:,2),in_situ_C.flow(:,1),in_situ_C.flow(:,2),'.-r')
% quiver(in_situ_B.pos(:,1),in_situ_B.pos(:,2),in_situ_B.flow(:,1),in_situ_B.flow(:,2),'.-r')
% xlim(lims_mlon); ylim(lims_mlat)
% % xlim([262.5,276.8]); ylim([64.0,68.7])
% pbaspect([range(x2),range(x3),1])

%
% [phi,mlon,mlat,E2_bg,E3_bg,v2_int,v3_int] = tools.replicate(in_situ_A,image,xg ...
%     ,flow_smoothing_window=4,show_plots=1,do_scale=false,do_rotate=false, ...
%     edge_method="contour",arc_definition="conductance",contour_values=[4.4,4.4] ...
%     ,harmonic_mask=[1,0,1]*50e3,fit_harmonic=true);
% [phi,mlon,mlat,E2_bg,E3_bg] = jules.tools.replicate(in_situ_B,image,xg, ...
%     flow_smoothing_window = 16, ...
%     show_plots = true, ...
%     save_plots = false, ...
%     direc = 'plots', ...
%     add_phi_background = false, ...
%     fit_harmonic = false, ...
%     num_replications = 512, ...
%     arc_definition = "conductance", ...
%     edge_method = "contour", ...
%     do_rotate = true, ...
%     do_scale = false, ...
%     contour_values = [13,13], ...
%     harmonic_mask = [1,0,1]*50e3, ...
%     flow_bg = [-120,-383] ...
%     );
