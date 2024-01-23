path = 'data\20230319.001_lp_1min-fitcal-vvels_lat.h5';
time = mean(h5read(path,'/Time/UnixTime'));
tid = 76;
aid = 13;
fprintf('Time = %s\n',datetime(time(tid),'ConvertFrom','posixtime'))
galt = h5read(path,'/VvelsGeoCoords/Altitude');
glon = h5read(path,'/VvelsGeoCoords/Longitude');
glat = h5read(path,'/VvelsGeoCoords/Latitude');
gv = h5read(path,'/VvelsGeoCoords/Velocity');
fprintf('Altitude = %i m\n',galt(1,aid));

glon = squeeze(glon(:,aid))';
glat = squeeze(glat(:,aid))';
gv_east = squeeze(gv(1,:,aid,tid))/1e3; % see Laundal & Richmond (2016), Eq. 57-59
gv_north = squeeze(gv(2,:,aid,tid))/1e3;
gv_up = squeeze(gv(3,:,aid,tid))/1e3;

figure
hold on
quiver(glon,glat,gv_east,gv_north,0,'.-')
scatter(-147,66)
xlim([-160,-135]); ylim([60,70])
grid