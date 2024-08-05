direc = '\\dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\swop_0210_AC_02';
direc_split = strsplit(direc,'_');
md = direc_split{end-2};
sats = direc_split{end-1};
cfg = gemini3d.read.config(direc);
% xg = gemini3d.read.grid(direc);

% model space
RE = 6371;
model_lon = double(xg.glon-7);
model_lat = double(xg.glat);
model_rad = double(xg.r);
model_az = deg2rad(model_lon);
model_el = deg2rad(model_lat);
model_ra = model_rad / 1e3;
[model_x,model_y,model_z] = sph2cart(model_az,model_el,model_ra);

model_mlat = squeeze(90-double(xg.theta(1,:,:))*180/pi);
model_mlon = squeeze(double(xg.phi(1,:,:)))*180/pi;
flon = griddedInterpolant(model_mlon,model_mlat,squeeze(model_lon(1,:,:)));
flat = griddedInterpolant(model_mlon,model_mlat,squeeze(model_lat(1,:,:)));

% Alaska
shp_usa = shaperead(fullfile('data','map_data','gadm41_USA_1.shp'));
shp_can = shaperead(fullfile('data','map_data','gadm41_CAN_1.shp'));

lim.lon = [-175,-100]+360;
lim.lat = [55,80];

map_lon_usa = [shp_usa.X]'+360;
map_lat_usa = [shp_usa.Y]';
map_lon_usa(map_lon_usa > lim.lon(2)) = NaN;
map_lon_usa(map_lon_usa < lim.lon(1)) = NaN;
map_lat_usa(map_lat_usa > lim.lat(2)) = NaN;
map_lat_usa(map_lat_usa < lim.lat(1)) = NaN;
map_theta_usa = (90-map_lat_usa)*pi/180;
map_phi_usa = map_lon_usa*pi/180;
map_x_usa = RE*cos(map_phi_usa).*sin(map_theta_usa);
map_y_usa = RE*sin(map_phi_usa).*sin(map_theta_usa);
map_z_usa = RE*cos(map_theta_usa);

map_lon_can = [shp_can.X]'+360;
map_lat_can = [shp_can.Y]';
map_lon_can(map_lon_can > lim.lon(2)) = NaN;
map_lon_can(map_lon_can < lim.lon(1)) = NaN;
map_lat_can(map_lat_can > lim.lat(2)) = NaN;
map_lat_can(map_lat_can < lim.lat(1)) = NaN;
map_theta_can = (90-map_lat_can)*pi/180;
map_phi_can = map_lon_can*pi/180;
map_x_can = RE*cos(map_phi_can).*sin(map_theta_can);
map_y_can = RE*sin(map_phi_can).*sin(map_theta_can);
map_z_can = RE*cos(map_theta_can);

% earth
[earth_x,earth_y,earth_z] = sphere(50);
earth_x = RE*earth_x;
earth_y = RE*earth_y;
earth_z = RE*earth_z;

%% trajectory
swarm_path = fullfile('data','swop','swarm_data');
swarm_fn = {dir(fullfile(swarm_path,['*2023',md,'*.h5'])).name};
for fn = swarm_fn
    f = fn{1};
    sat = f(12);
    swarm_ids.(sat) = double(h5read(fullfile(swarm_path,f),'/ids'));
    swarm_time.(sat) = datetime(h5read(fullfile(swarm_path,f),'/Timestamp'),'convertfrom','posixtime');
    swarm_lon.(sat) = h5read(fullfile(swarm_path,f),'/Longitude');
    swarm_lat.(sat) = h5read(fullfile(swarm_path,f),'/Latitude');
    swarm_vie.(sat) = h5read(fullfile(swarm_path,f),'/ViE') / 1e3;
    swarm_vin.(sat) = h5read(fullfile(swarm_path,f),'/ViN') / 1e3;
    swarm_rad.(sat) = h5read(fullfile(swarm_path,f),'/Radius') / 1e3;
end

for s = sats
    ids = swarm_ids.(s);
    id0 = ids(1);
    id1 = ids(2);
    lon = swarm_lon.(s);
    disp(mean(lon(id0:id1))+360)
    lat = swarm_lat.(s);
    rad = swarm_rad.(s);
    az = deg2rad(lon(id0:id1));
    el = deg2rad(lat(id0:id1));
    ra = rad(id0:id1);
    [traj_x.(s),traj_y.(s),traj_z.(s)] = sph2cart(az,el,ra);
    [traj_shd_x.(s),traj_shd_y.(s),traj_shd_z.(s)] = sph2cart(az,el,RE);
end

% points
min_x = model_x(1,end,end);
min_y = model_y(1,end,end);
min_z = model_z(1,end,end);
max_x = model_x(end,end,end);
max_y = model_y(end,end,end);
max_z = model_z(end,end,end);

[fair_x,fair_y,fair_z] = sph2cart(deg2rad(-147.7200+360),deg2rad(64.8401),RE+9);
[anch_x,anch_y,anch_z] = sph2cart(deg2rad(-149.8997+360),deg2rad(61.2176),RE+9);

% plot 3d context
close all
fts = 18*2;
ftn = 'Arial';
lw = 2;

reset(0)
set(0,'defaultFigurePaperUnits','inches')
jules.tools.setall(0,'FontName',ftn)
jules.tools.setall(0,'FontSize',fts)
jules.tools.setall(0,'Multiplier',1)
set(0,'defaultAxesFontSizeMode','manual')
set(0,'defaultSurfaceEdgeColor','flat')
set(0,'defaultLineLineWidth',lw)
set(0,'defaultScatterLineWidth',lw)
set(0,'defaultQuiverLineWidth',lw)

vv = [-35,0];
cz = 5.5;
cp = [-0.4,0.7];
ar = [1,1,1];

figure
set(gcf,'PaperPosition',[0,0,10,10])

hold on
axis off
cc = [30,92,141]/255;
plot3(map_x_usa,map_y_usa,map_z_usa,'Color',[1,1,1]*0.4)
plot3(map_x_can,map_y_can,map_z_can,'Color',[1,1,1]*0.4)
model(model_x,model_y,model_z,[0,0.7,0])
text(max_x-1100,max_y,max_z+30,'Model space','Color',[0,0.7,0])
scatter3(min_x,min_y,min_z,200,'filled','b')
text(min_x+100,min_y,min_z ...
    ,[num2str(round(min(model_rad(:)/1e3))-RE),' km'],'Color','b')
scatter3(max_x,max_y,max_z,200,'filled','b')
text(max_x+100,max_y,max_z ...
    ,[num2str(round(max(model_rad(:)/1e3))-RE),' km'],'Color','b')
for s = sats
    plot3(traj_x.(s),traj_y.(s),traj_z.(s),'Color','r')
    plot3(traj_shd_x.(s),traj_shd_y.(s),traj_shd_z.(s),'--','Color',[1,0.5,0.5])
end
text(min_x-730,min_y,min_z-450,'Swarm A','Color','r')
text(min_x-1300,min_y,min_z-490,'Swarm C','Color','r')
scatter3(fair_x,fair_y,fair_z,200,'filled','MarkerFaceColor',cc)
text(fair_x-450,fair_y,fair_z+40,'Fairbanks','Color',cc,'FontSize',fts*0.8)
scatter3(anch_x,anch_y,anch_z,200,'filled','MarkerFaceColor',cc)
text(anch_x-530,anch_y,anch_z+40,'Anchorage','Color',cc,'FontSize',fts*0.8)
view(vv)
camzoom(cz)
camdolly(cp(1),cp(2),0)
pbaspect(ar)

saveas(gcf,fullfile('plots','swop','context_3d_swop.png'))
close all

%% load data maps
map = h5read(fullfile(direc,'inputs','particles', ...
    [char(gemini3d.datelab(cfg.times(end))),'.h5']),'/Qp');
mlon = h5read(fullfile(direc,'inputs','fields','simgrid.h5'),'/mlon');
mlat = h5read(fullfile(direc,'inputs','fields','simgrid.h5'),'/mlat');
[image_mlon,image_mlat] = ndgrid(mlon,mlat);
image_lon = flon(image_mlon,image_mlat);
image_lat = flat(image_mlon,image_mlat);

% plot 2d context
close all

colorcet = @jules.tools.colorcet;

figure(2)
set(gcf,'PaperPosition',[0,0,10,6])
sc = 0.2;
clm = 'L5';
lim.x = [min(image_lon(:)),max(image_lon(:))];
lim.y = [min(image_lat(:)),max(image_lat(:))];

hold on
pcolor(image_lon,image_lat,map)
plot(squeeze(model_lon(1,:,1)),squeeze(model_lat(1,:,1)),'Color',[0,0.7,0])
plot(squeeze(model_lon(1,:,end)),squeeze(model_lat(1,:,end)),'Color',[0,0.7,0])
plot(squeeze(model_lon(1,1,:)),squeeze(model_lat(1,1,:)),'Color',[0,0.7,0])
plot(squeeze(model_lon(1,end,:)),squeeze(model_lat(1,end,:)),'Color',[0,0.7,0])
% plot(traj_lon,traj_lat,'r','LineWidth',lw*1.5)
for s = sats
    ids = swarm_ids.(s);
    id0 = ids(1);
    id1 = ids(2);
    lon = swarm_lon.(s);
    lat = swarm_lat.(s);
    lon = lon(id0:id1);
    lat = lat(id0:id1);
    ve = swarm_vie.(s)';
    vn = swarm_vin.(s)';
    ids = 600:length(ve);
    quiver(lon(ids)+360,lat(ids),ve(ids)*sc,vn(ids)*sc,0,'.-r','LineWidth',lw/2)
end
pbaspect([1.618,1,1])
xlim(lim.x); ylim(lim.y)
xlabel('Geo. longitude (°)'); ylabel('Geo. latitude (°)')
grid on
colormap(gca,colorcet(clm))
% clb = colorbar;
% clb.Label.String = 'e^- Flux (mW/m^2)';
camdolly(-0.2,-0.1,0)
% camzoom(1)

saveas(gcf,fullfile('plots','swop','context_2d_swop.png'))
close all

%% combine images
% fileA = fullfile('figures','context_3d_arcs.png');
% fileB = fullfile('figures','context_3d_ig.png');
% fileC = fullfile('figures','context_2d_arcs.png');
% fileD = fullfile('figures','context_2d_ig.png');
% if all([isfile(fileA),isfile(fileB),isfile(fileC),isfile(fileD)])
%     imA = imread(fileA);
%     imB = imread(fileB);
%     imC = imread(fileC);
%     imD = imread(fileD);
%     im = cat(1,cat(2,imA,imB),cat(2,imC,imD));
%     im = insertText(im,[0,0] ...
%         ,'A','FontSize',40,'Font','Arial','TextBoxColor','w');
%     im = insertText(im,[size(imA,2),0] ...
%         ,'B','FontSize',40,'Font','Arial','TextBoxColor','w');
%     im = insertText(im,[0,size(imA,1)] ...
%         ,'C','FontSize',40,'Font','Arial','TextBoxColor','w');
%     im = insertText(im,[size(imA,2),size(imA,1)] ...
%         ,'D','FontSize',40,'Font','Arial','TextBoxColor','w');
%     imwrite(im,fullfile('figures','01_context.png'))
% end

%% model space edges
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