ig = true;
%#ok<*UNRCH>

if ig
    direc = '//dartfs-hpc/rc/lab/L/LynchK/public_html/Gemini3D/isinglass_74';
else
    direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\aurora_longsharc_bent_01';
end
xg = gemini3d.read.grid(direc);

% model space
RE = 6371;
model_lon = double(xg.glon);
model_lat = double(xg.glat);
model_phi = model_lon*pi/180;
model_theta = (90-model_lat)*pi/180;
model_alt = double(xg.alt)/1e3 + RE;
model_x = model_alt.*cos(model_phi).*sin(model_theta);
model_y = model_alt.*sin(model_phi).*sin(model_theta);
model_z = model_alt.*cos(model_theta);

model_mlat = squeeze(90-double(xg.theta(1,:,:))*180/pi);
model_mlon = squeeze(double(xg.phi(1,:,:)))*180/pi;
flon = griddedInterpolant(model_mlon,model_mlat,squeeze(model_lon(1,:,:)));
flat = griddedInterpolant(model_mlon,model_mlat,squeeze(model_lat(1,:,:)));

% Alaska
shp = shaperead('data/map_data/gadm41_USA_2.shp');

lim.lon = [-170,-141]+360;
lim.lat = [55,72];
map_lon = [shp.X]'+360;
map_lat = [shp.Y]';
map_lon(map_lon > lim.lon(2)) = NaN;
map_lon(map_lon < lim.lon(1)) = NaN;
map_lat(map_lat > lim.lat(2)) = NaN;
map_lat(map_lat < lim.lat(1)) = NaN;

map_theta = (90-map_lat)*pi/180;
map_phi = map_lon*pi/180;
map_x = RE*cos(map_phi).*sin(map_theta);
map_y = RE*sin(map_phi).*sin(map_theta);
map_z = RE*cos(map_theta);

% earth
[earth_x,earth_y,earth_z] = sphere(50);
earth_x = RE*earth_x;
earth_y = RE*earth_y;
earth_z = RE*earth_z;

%% trajectory
if ig
    load('data\flow_data_ig.mat','mlon','mlat','v_geog_east','v_geog_north')
    traj_lon = flon(mlon,mlat);
    traj_lat = flat(mlon,mlat);
    data = strsplit(fileread('data\GPSecef.txt'),'\n');
    traj_x = nan(1,length(data)-22);
    traj_y = traj_x;
    traj_z = traj_x;
    j=0;
    for i = 22:length(data)-1
        datum = strsplit(data{i});
        traj_x(i-21) = str2double(cell2mat(datum(3)))/1e3;
        traj_y(i-21) = str2double(cell2mat(datum(4)))/1e3;
        traj_z(i-21) = str2double(cell2mat(datum(5)))/1e3;
        j=j+1;
    end
    [az,el,~] = cart2sph(traj_x,traj_y,traj_z);
    [traj_shd_x,traj_shd_y,traj_shd_z] = sph2cart(az,el,RE);
    traj_vlon = smoothdata(v_geog_east,"gaussian",16);
    traj_vlat = smoothdata(v_geog_north,"gaussian",16);
else
    load(fullfile(direc,'arcs_orbit_20210211_track.mat'),'glonsat','glatsat','altsat')
    az = deg2rad(glonsat(1:100:end,4:4:end));
    el = deg2rad(glatsat(1:100:end,4:4:end));
    ra = altsat(1:100:end,4:4:end)/1e3 + RE;
    [traj_x,traj_y,traj_z] = sph2cart(az,el,ra);
    [traj_shd_x,traj_shd_y,traj_shd_z] = sph2cart(az,el,RE);
    traj_lon = rad2deg(az)+360;
    traj_lat = rad2deg(el);
end

%% points
min_x = model_x(1,end,end);
min_y = model_y(1,end,end);
min_z = model_z(1,end,end);
max_x = model_x(end,end,end);
max_y = model_y(end,end,end);
max_z = model_z(end,end,end);

[fair_x,fair_y,fair_z] = sph2cart(deg2rad(-147.7200+360),deg2rad(64.8401),RE+10);
[anch_x,anch_y,anch_z] = sph2cart(deg2rad(-149.8997+360),deg2rad(61.2176),RE+10);

% plot
close all
fts = 20;
ftn = 'Arial';
lw = 1.4;

reset(0)
set(0,'defaultFigurePaperUnits','inches')
tools.setall(0,'FontName',ftn)
tools.setall(0,'FontSize',fts)
tools.setall(0,'Multiplier',1)
set(0,'defaultAxesFontSizeMode','manual')
set(0,'defaultSurfaceEdgeColor','flat')
set(0,'defaultLineLineWidth',lw)
set(0,'defaultScatterLineWidth',lw)
set(0,'defaultQuiverLineWidth',lw)

if ig
    vv = [-30,10];
    cz = 1.5;
    cp = [0,0];
    ar = [1,1,1];
    A = [-1000,0,150];
    B = [-350,150,-550];
    C = [0,-10,0,-50,-500,80];
    trajectory = 'Isinglass trajectory';
else
    vv = [-45,10];
    cz = 1.8;
    cp = [0.5,0];
    ar = [1,1.3,1];
    A = [-1000,0,400];
    B = [-1900,-800,0];
    C = [0,min_x-traj_x(1),0,min_y-traj_y(1),0,min_z-traj_z(1)];
    trajectory = 'ARCS trajectories';
end

figure
set(gcf,'PaperPosition',[0,0,13.2,13.2]/2)

hold on
axis off
cc = [30,92,141]/255;
% surf(earth_x,earth_y,earth_z)
plot3(map_x,map_y,map_z,'Color',[1,1,1]*0.4)
model(model_x,model_y,model_z,[0,0.7,0])
text(max_x+A(1),max_y+A(2),max_z+A(3),'Model space','Color',[0,0.7,0])
scatter3(min_x,min_y,min_z,'filled','b')
text(min_x+100,min_y,min_z,[num2str(round(min(model_alt(:)))-RE),' km'],'Color','b')
scatter3(max_x,max_y,max_z,'filled','b')
text(max_x+100,max_y,max_z,[num2str(round(max(model_alt(:)))-RE),' km'],'Color','b')
plot3(traj_x,traj_y,traj_z,'Color','r')
plot3(traj_shd_x,traj_shd_y,traj_shd_z,'--','Color',[1,0.5,0.5])
text(min_x+B(1),min_y+B(2),min_z+B(3),trajectory,'Color','r')
plot3([min_x+C(1),traj_x(1)+C(2)] ...
    ,[min_y+C(3),traj_y(1)+C(4)] ...
    ,[min_z+C(5),traj_z(1)+C(6)],'Color','r','LineWidth',0.5*lw)
if ig
    scatter3(fair_x,fair_y,fair_z,100,'filled','MarkerFaceColor',cc)
    text(fair_x-390,fair_y,fair_z+70,'Fairbanks','Color',cc,'FontSize',fts*0.8)
    scatter3(anch_x,anch_y,anch_z,100,'filled','MarkerFaceColor',cc)
    text(anch_x-450,anch_y,anch_z+70,'Anchorage','Color',cc,'FontSize',fts*0.8)
end
view(vv)
camzoom(cz)
campan(cp(1),cp(2))
pbaspect(ar)

if ig
    saveas(gcf,fullfile('plots','paper0','context_ig.png'))
else
    saveas(gcf,fullfile('plots','paper0','context_arcs.png'))
end
close all

%%
if ig
    load('data\particles.mat','Qit','E0it','mlat','mlon')
    map = Qit(:,:,4);
else
    map = h5read(fullfile(direc,'inputs','fields','20150201_36000.000000.h5'),'/Vmaxx1it')*1e6;
    mlon = h5read(fullfile(direc,'inputs','fields','simgrid.h5'),'/mlon');
    mlat = h5read(fullfile(direc,'inputs','fields','simgrid.h5'),'/mlat');
end
[image_mlon,image_mlat] = ndgrid(mlon,mlat);
image_lon = flon(image_mlon,image_mlat);
image_lat = flat(image_mlon,image_mlat);

%%
close all

figure(2)
set(gcf,'PaperPosition',[0,0,13.2/2,5])
sc = 3e-4;

if ig
    ids = 1:2:length(traj_vlat);
    xlims = [210.6,216.7];
    ylims = [65.2,67.3];
    clims = [1,49];
    clm = 'kbgyw';
    clb_lbl = 'Q (mW/m^2)';
else
    xlims = [165,240];
    ylims = [55,77];
    clims = [-1,1];
    clm = 'D1A';
    clb_lbl = 'j_{||} (uA/m^2)';
end

hold on
pcolor(image_lon,image_lat,map)
colormap(gca,colorcet(clm))
% plot(map_lon,map_lat,'k')
plot(squeeze(model_lon(1,:,1)),squeeze(model_lat(1,:,1)),'Color',[0,0.7,0])
plot(squeeze(model_lon(1,:,end)),squeeze(model_lat(1,:,end)),'Color',[0,0.7,0])
plot(squeeze(model_lon(1,1,:)),squeeze(model_lat(1,1,:)),'Color',[0,0.7,0])
plot(squeeze(model_lon(1,end,:)),squeeze(model_lat(1,end,:)),'Color',[0,0.7,0])
plot(traj_lon,traj_lat,'r','LineWidth',lw*1.5)
if ig
    quiver(traj_lon(ids),traj_lat(ids),traj_vlon(ids)*sc,traj_vlat(ids)*sc+0.05,0,'.-r','LineWidth',lw/2)
end
pbaspect([1,1,1])
xlim(xlims); ylim(ylims); clim(clims)
xlabel('G. longitude (°)'); ylabel('G. latitude (°)')
grid on
clb = colorbar;
clb.Label.String = clb_lbl;

if ig
    saveas(gcf,fullfile('plots','paper0','context2_ig.png'))
else
    saveas(gcf,fullfile('plots','paper0','context2_arcs.png'))
end
close all

%%
function model(x,y,z,c)
lx1 = size(x,1);
lx2 = size(x,2);
lx3 = size(x,3);
for i1 = [1,lx1]
    for i2 = [1,lx2]
        xx = squeeze(x(i1,i2,:));
        yy = squeeze(y(i1,i2,:));
        zz = squeeze(z(i1,i2,:));
        plot3(xx,yy,zz,'Color',c)
    end
end
for i1 = [1,lx1]
    for i3 = [1,lx3]
        xx = squeeze(x(i1,:,i3));
        yy = squeeze(y(i1,:,i3));
        zz = squeeze(z(i1,:,i3));
        plot3(xx,yy,zz,'Color',c)
    end
end
for i2 = [1,lx2]
    for i3 = [1,lx3]
        xx = squeeze(x(:,i2,i3));
        yy = squeeze(y(:,i2,i3));
        zz = squeeze(z(:,i2,i3));
        plot3(xx,yy,zz,'Color',c)
    end
end
end