%% read config + grid data
direc = '../../public_html/Gemini3D/isinglass_73/';
cfg = gemini3d.read.config(direc);
xg = gemini3d.read.grid(direc);
dat = gemini3d.read.frame(direc,'time',cfg.times(end));

%% unpack grid data
x2 = xg.x2(3:end-2);
x3 = xg.x3(3:end-2);
[X2,X3] = ndgrid(x2,x3);
MLAT = 90-squeeze(xg.theta(end,:,:))*180/pi;
MLON = squeeze(xg.phi(end,:,:))*180/pi;
mlon_to_x2 = griddedInterpolant(MLON(:,1),xg.x2(3:end-2));
mlat_to_x3 = griddedInterpolant(MLAT(1,:),xg.x3(3:end-2));

% read isinglass data
fac = -squeeze(dat.J1(end,:,:))*1e6;
ffac = griddedInterpolant(X2,X3,fac,'linear');
% load(fullfile(direc,'ext','flow_data.mat'),'mlon','mlat','time')
load('data\flow_data.mat','mlon','mlat','time')
x2_isin = mlon_to_x2(mlon);
x3_isin = mlat_to_x3(mlat);
buf = 0.97;
ids = x3_isin > min(x3)*buf & x3_isin < max(x3)*buf;
x2_isin = x2_isin(ids);
x3_isin = x3_isin(ids);
time = time(ids);
fac_isin = ffac(x2_isin,x3_isin);

load('data\isinglass_fac_data.mat','jpartime','jpar')

%% plot
close all

j_lim = 120;
clm = colorcet('D1A');
co = colororder;

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

% figure(1)
% hold on
% % pcolor(X2/1e3,X3/1e3,smoothdata(fac,2,'gaussian',6))
% pcolor(X2/1e3,X3/1e3,fac)
% plot(x2_isin/1e3,x3_isin/1e3,'m')
% shading flat
% xlabel('East (km)')
% ylabel('North (km)')
% clb = colorbar;
% clb.Label.String = 'fac (\muA/m^2)';
% clim([-1,1]*j_lim*2)
% colormap(clm)

figure(2)
set(gcf,'PaperPosition',[0,0,13,4.1])

hold on
plot(jpartime,jpar,'k')
plot(time,fac_isin,':r')
plot(1.13*(time-time(17))+time(17),fac_isin,'r')
xlim([150,400])
ylim([-1,1]*j_lim)
xlabel('Flight time (s)')
ylabel('Field-alligned current (\muA/m^2)')
grid
legend('Isinglass in situ','Model','Model scaled by 13%')

% saveas(gcf,'plots\fac_comparison_IG_x1.1-scaled.png')
% saveas(gcf,'plots\fac_comparison_IG_12s-delay.png')
saveas(gcf,'plots\fac_comparison_IG.png')
% close all
