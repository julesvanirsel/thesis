clear
md = '0314';
sats = 'C';
ver = '02';
dasc_path = fullfile('data','swop','dasc_data');
swarm_path = fullfile('data','swop','swarm_data');
dasc_fn = dir(fullfile(dasc_path,['2023',md,'*.mat'])).name;
swarm_fn = {dir(fullfile(swarm_path,['*2023',md,'*.h5'])).name};
clr.A = 'r'; clr.B = 'b'; clr.C = 'c';

load(fullfile(dasc_path,dasc_fn),'data')
dasc_data = data;
clear('data')
dasc_mlon = dasc_data.maglon_110 + 360;
dasc_mlat = dasc_data.maglat_110;
dasc_sigp = dasc_data.SigP;
dasc_sigh = dasc_data.SigH;

for fn = swarm_fn
    f = fn{1};
    sat = f(12);
    swarm_ids.(sat) = double(h5read(fullfile(swarm_path,f),'/ids'));
    swarm_time.(sat) = datetime(h5read(fullfile(swarm_path,f),'/Timestamp'),'convertfrom','posixtime');
    swarm_mlon.(sat) = h5read(fullfile(swarm_path,f),'/MagneticLongitude');
    swarm_mlat.(sat) = h5read(fullfile(swarm_path,f),'/MagneticLatitude');
    swarm_vime.(sat) = h5read(fullfile(swarm_path,f),'/ViMagE') / 1e3;
    swarm_vimn.(sat) = h5read(fullfile(swarm_path,f),'/ViMagN') / 1e3;
    swarm_gdal.(sat) = h5read(fullfile(swarm_path,f),'/GeodeticAltitude') / 1e3;
end

% sats = cell2mat(fieldnames(swarm_mlon))';

% alex_path = fullfile('data','swop','flow_data_alex');
% swarm_fn_alex = {dir(fullfile(alex_path,['*2023',mds,'*.mat'])).name};
% dlon = -0.4*0;
% dlat = 0.0394*0;
% ltr = 'A';
% for fn = swarm_fn_alex
%     f = fn{1};
%     load(fullfile(alex_path,f))
%     alex_data = eval(f(1:end-17));
%     clear(f(1:end-10))
%     alex_time_tmp = string(f(20:27)) + pad(string(alex_data.UTtime),10,'left','0');
%     alex_time.(ltr) = datetime(alex_time_tmp,'inputformat','yyyyMMddHHmmss.SSS');
%     alex_mlon.(ltr) = alex_data.Mag_Lon + 360 + dlon;
%     alex_mlat.(ltr) = alex_data.Mag_Lat + dlat;
%     alex_vime.(ltr) = alex_data.Flow_MagEast / 1e3;
%     alex_vimn.(ltr) = alex_data.Flow_MagNorth / 1e3;
%     alex_gdal.(ltr) = alex_data.Geodetic_Altitude_km;
%     ltr = char(ltr+1);
% end

xlims = [min(dasc_mlon(:)),max(dasc_mlon(:))];
ylims = [min(dasc_mlat(:)),max(dasc_mlat(:))];

lims_all.Dnull = [[min(dasc_mlon(:)), max(dasc_mlon(:))]; ...
    [min(dasc_mlat(:)), max(dasc_mlat(:))]];
lims_all.D0210 = [[268.1, 275.8]; [65.00, 67.05]];
lims_all.D0212 = lims_all.Dnull;
lims_all.D0214 = lims_all.Dnull;
lims_all.D0304 = [[259.3, 266.5]; [64.4, 66.4]];
lims_all.D0314 = [[262.5, 276.8]; [62.7, 69]];
lims_all.D0319 = [[264.9, 277.1]; [65.70, 68.70]];
ctrs_all.D0210 = 16;
ctrs_all.D0212 = 11;
ctrs_all.D0214 = 11;
ctrs_all.D0304 = 5;
ctrs_all.D0314 = 5;
ctrs_all.D0319 = 16;
arc_all.D0210 = dasc_sigh;
arc_all.D0212 = dasc_sigh;
arc_all.D0214 = dasc_sigh;
arc_all.D0304 = dasc_sigh;
arc_all.D0314 = dasc_sigh;
arc_all.D0319 = dasc_sigh;
lims = lims_all.(['D',md]);
ctrs = ctrs_all.(['D',md]);
arc = smoothdata2(arc_all.(['D',md]),"gaussian",30);

clims = [quantile(arc(:),0.01), quantile(arc(:),0.99)];

for sat = sats
    ok.(sat) = (swarm_mlat.(sat) > lims(2,1)) & (swarm_mlat.(sat) < lims(2,2));
end

% idsA = 560:747;
% idsC = 550:737;
% vbg = [108.6123, 468.4083]/1e3;
% t = 30;
% vbg = [sin(deg2rad(t)), cos(deg2rad(t))];
% vbg = [700,-900]/1e3;

close all
farc = griddedInterpolant(dasc_mlon,dasc_mlat,arc);
for sat = sats
    arc_cut = farc(swarm_mlon.(sat), swarm_mlat.(sat));
    [~,vbg_id.(sat)] = max(arc_cut);
end

close all
figure('Position',[10 50 1500 800])
hold on
title(sprintf('md = %s   sats = %s',md,sats))
pcolor(dasc_mlon,dasc_mlat,arc)
% contour(dasc_mlon,dasc_mlat,arc,'w')
ctr = contour(dasc_mlon,dasc_mlat,arc,[1,1]*ctrs,'g');
% [x0,y0] = jules.tools.unpack_contours(ctr,rank=1);
% [x1,y1] = jules.tools.unpack_contours(ctr,rank=-1);
% plot(x0,y0,'--r')
% plot(x1,y1,'--r')
for sat = sats
    mlon_tmp = swarm_mlon.(sat);
    mlat_tmp = swarm_mlat.(sat);
    vime_tmp = swarm_vime.(sat);
    vimn_tmp = swarm_vimn.(sat);
    % vbg = [vime_tmp(vbg_id.(sat)), vimn_tmp(vbg_id.(sat))];
    vbg = [0,0];
    oks = ok.(sat);
    c = clr.(sat);
    quiver(mlon_tmp(oks),mlat_tmp(oks),vime_tmp(oks)-vbg(1),vimn_tmp(oks)-vbg(2),0,['.-',c])
end
shading flat
colorbar;
colormap('pink')
xlim(lims(1,:)); ylim(lims(2,:)); clim(clims)
pbaspect([range(lims(1,:))*0.4, range(lims(2,:)), 1])

%% from gemini3d.grid.cartesian.m
Re=6370e3;

theta = deg2rad(90-lims(2,:));
thetactr = mean(theta);
gamma2 = thetactr - theta;
y = Re*gamma2;
ydist = range(y);
ydist = ceil(ydist/1e4)*1e4 + 10e3;

phi = deg2rad(lims(1,:));
phictr = mean(phi);
gamma1 = phi - phictr;
x = gamma1*sin(thetactr)*Re;
xdist = range(x);
xdist = ceil(xdist/1e4)*1e4 + 10e3;

fprintf('%s\n',pad(md,80,'both','-'))
jules.tools.grid_params([320,448],[xdist,ydist]/1e3,do_plot=true);

[glat,glon] = gemini3d.geomag2geog(thetactr,phictr);
fprintf('glat = %.3f\n',glat)
fprintf('glon = %.2f\n',glon)
fprintf('xdist = %.0fe3\n',xdist/1e3)
fprintf('ydist = %.0fe3\n',ydist/1e3)

%%
direc = ['/lab/L/LynchK/public_html/Gemini3D/swop_',md,'_',sats(1),'_',ver];
cfg = gemini3d.read.config(direc);
xg = gemini3d.grid.cartesian(cfg);

%% check grid
mlon = squeeze(rad2deg(xg.phi(1,:,1)));
mlat = squeeze(90-rad2deg(xg.theta(1,1,:)))';

close all
figure('Position',[10 50 1500 800])
hold on
pcolor(dasc_mlon,dasc_mlat,arc)
contour(dasc_mlon,dasc_mlat,arc,[1,1]*ctrs,'g')
for sat = sats
    mlon_tmp = swarm_mlon.(sat);
    mlat_tmp = swarm_mlat.(sat);
    vime_tmp = swarm_vime.(sat);
    vimn_tmp = swarm_vimn.(sat);
    oks = ok.(sat);
    c = clr.(sat);
    quiver(mlon_tmp(oks),mlat_tmp(oks),vime_tmp(oks),vimn_tmp(oks),0,['.-',c])
end
shading flat
colorbar;
colormap('pink')
xlim([min(mlon),max(mlon)]); ylim([min(mlat),max(mlat)]); clim(clims)
pbaspect([range(lims(1,:))*0.4, range(lims(2,:)), 1])
