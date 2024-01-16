%%
close all
direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\swop_01';
cfg = gemini3d.read.config(direc);
xg = gemini3d.grid.cartesian(cfg);

%%
mlon = squeeze(xg.phi(1,:,1)*180/pi);
mlat = squeeze(90-xg.theta(1,1,:)*180/pi);
x2 = double(xg.x2(3:end-2));
x3 = double(xg.x3(3:end-2));

mlon_to_x2 = griddedInterpolant(mlon,x2);
mlat_to_x3 = griddedInterpolant(mlat,x3);

lims_mlon = [min(mlon),max(mlon)];
lims_mlat = [min(mlat),max(mlat)];

%%
% load('data\replicate_data_isinglass.mat','image','in_situ')
% image_IG = image;
% in_situ_IG = in_situ;
% clear('image','in_situ')

load('data\0314_precip_map.mat')
load('data\SW_EXPT_EFIA_TCT02_20230314T064200.mat')
load('data\SW_EXPT_EFIC_TCT02_20230314T064200.mat')

% load('data\0212_precip_map_lowres.mat')
% load('data\SW_EXPT_EFIC_TCT02_20230212T101600.mat')

% lims_mlon = [268.84, 276.26];
% lims_mlat = [63.65, 67.70];
event = dt_20230314T0649;

lims_mlon = [min(event.maglon_110(:))+360,max(event.maglon_110(:))+360];
lims_mlat = [min(event.maglat_110(:)),max(event.maglat_110(:))];

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
% clear('event')

% figure
% tiledlayout(2,2)
% nexttile
% pcolor(image.pos(:,:,1)); colorbar
% nexttile
% pcolor(image.pos(:,:,2)); colorbar
% nexttile
% pcolor(image_IG.pos(:,:,1)); colorbar
% nexttile
% pcolor(image_IG.pos(:,:,2)); colorbar

%%
SW_EXPT_A = SW_EXPT_EFIA_TCT02_20230314T064200;
SW_EXPT_C = SW_EXPT_EFIC_TCT02_20230314T064200;
% SW_EXPT_C = SW_EXPT_EFIC_TCT02_20230212T101600;

in_situ_A.pos(:,1) = SW_EXPT_A.Mag_Lon + 360;
in_situ_A.pos(:,2) = SW_EXPT_A.Mag_Lat;
in_situ_A.flow(:,1) = SW_EXPT_A.Flow_MagEast;
in_situ_A.flow(:,2) = SW_EXPT_A.Flow_MagNorth;
% clear('SW_EXPT_A')

in_situ_C.pos(:,1) = SW_EXPT_C.Mag_Lon+360;
in_situ_C.pos(:,2) = SW_EXPT_C.Mag_Lat;
in_situ_C.flow(:,1) = SW_EXPT_C.Flow_MagEast;
in_situ_C.flow(:,2) = SW_EXPT_C.Flow_MagNorth;
% clear('SW_EXPT_C')

%%
% save('data\replicate_data_swop_01.mat','image','in_situ_C')
save('data\replicate_data_swop_02.mat','image','in_situ_A','in_situ_C')

%%
close all

figure
hold on
pcolor(image.pos(:,:,1),image.pos(:,:,2),image.flux); colorbar
quiver(in_situ_A.pos(:,1),in_situ_A.pos(:,2),in_situ_A.flow(:,1),in_situ_A.flow(:,2),'.-r')
quiver(in_situ_C.pos(:,1),in_situ_C.pos(:,2),in_situ_C.flow(:,1),in_situ_C.flow(:,2),'.-r')
xlim(lims_mlon); ylim(lims_mlat)
pbaspect([range(x2),range(x3),1])

figure
hold on
pcolor(image.pos(:,:,1),image.pos(:,:,2),image.energy); colorbar
quiver(in_situ_A.pos(:,1),in_situ_A.pos(:,2),in_situ_A.flow(:,1),in_situ_A.flow(:,2),'.-r')
quiver(in_situ_C.pos(:,1),in_situ_C.pos(:,2),in_situ_C.flow(:,1),in_situ_C.flow(:,2),'.-r')
xlim(lims_mlon); ylim(lims_mlat)
pbaspect([range(x2),range(x3),1])

%%
figure
set(gcf,'PaperPosition',[0,0,6.5,5.5])
set(gca,'Fontsize',20)
sc = 0.001;
sc2 = 0.0005;
ids = 1:3:length(in_situ_A.flow);

hold on
pcolor(image.pos(:,:,1),image.pos(:,:,2),image.flux); shading flat; colormap(colorcet('L18'))
contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,10)
quiver(in_situ_A.pos(ids,1),in_situ_A.pos(ids,2),in_situ_A.flow(ids,1)*sc,in_situ_A.flow(ids,2)*sc,0,'.-b')
quiver(in_situ_C.pos(ids,1),in_situ_C.pos(ids,2),in_situ_C.flow(ids,1)*sc2,in_situ_C.flow(ids,2)*sc2,0,'.-b')
xlim(lims_mlon); ylim(lims_mlat)
xlabel('M. Longitude'); ylabel('M. Latitude')
clb = colorbar;
clb.Label.String = 'Q (mW/m^2)  and  \Sigma_P (S)';
clim([0,10])
pbaspect([1,1,1])

saveas(gcf,'plots/swop_dasc.png')

%%
disp(any(isnan(image.flux(:))))
disp(any(isnan(image.energy(:))))
disp(any(isnan(image.pedersen(:))))

%%
figure
hold on
contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,20); colorbar; shading flat
contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,[1,1]*4.5,'m'); colorbar; shading flat
contour(image.pos(:,:,1),image.pos(:,:,2),image.pedersen,[1,1]*4.5,'m'); colorbar; shading flat
% quiver(in_situ_A.pos(:,1),in_situ_A.pos(:,2),in_situ_A.flow(:,1),in_situ_A.flow(:,2),'.-r')
quiver(in_situ_C.pos(:,1),in_situ_C.pos(:,2),in_situ_C.flow(:,1),in_situ_C.flow(:,2),'.-r')
xlim(lims_mlon); ylim(lims_mlat)
% xlim([262.5,276.8]); ylim([64.0,68.7])
pbaspect([range(x2),range(x3),1])

%%
% [phi,mlon,mlat,E2_bg,E3_bg,v2_int,v3_int] = tools.replicate(in_situ_A,image,xg ...
%     ,flow_smoothing_window=4,show_plots=1,do_scale=false,do_rotate=false, ...
%     edge_method="contour",arc_definition="conductance",contour_values=[4.4,4.4] ...
%     ,harmonic_mask=[1,0,1]*50e3,fit_harmonic=true);
[phi,mlon,mlat,E2_bg,E3_bg] = tools.replicate(in_situ_A,image,xg, ...
    flow_smoothing_window = 16, ...
    show_plots = true, ...
    save_plots = false, ...
    direc = 'plots', ...
    add_phi_background = false, ...
    fit_harmonic = false, ...
    num_replications = 512, ...
    arc_definition = "conductance", ...
    edge_method = "contour", ...
    do_rotate = true, ...
    do_scale = false, ...
    contour_values = [4.4,4.4], ...
    harmonic_mask = [1,0,1]*50e3, ...
    flow_bg = [-120,-383] ...
    );
