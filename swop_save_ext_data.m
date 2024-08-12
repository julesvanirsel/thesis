clear
direc = '\\dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\swop_0314_AC_02';
dasc_path = fullfile('data','swop','dasc_data');
swarm_path = fullfile('data','swop','swarm_data');
mkdir(fullfile(direc,'ext'))
md = direc(56:59);

cfg = gemini3d.read.config(direc);
xg = gemini3d.grid.cartesian(cfg);
gem_mlon = squeeze(rad2deg(xg.phi(1,:,1)));
gem_mlat = squeeze(90-rad2deg(xg.theta(1,1,:)))';

h5make = @jules.tools.h5make;

%% imagery data
dasc_fn = dir(fullfile(dasc_path,['2023*',md,'*.mat'])).name;
load(fullfile(dasc_path,dasc_fn),'data')
dasc = data;
clear('data')

[~,id00] = min(abs(gem_mlon(1)-dasc.maglon_110(:,1)-360));
[~,id01] = min(abs(gem_mlon(end)-dasc.maglon_110(:,1)-360));
[~,id10] = min(abs(gem_mlat(1)-dasc.maglat_110(1,:)));
[~,id11] = min(abs(gem_mlat(end)-dasc.maglat_110(1,:)));

% sp = 18; %0210
sp = 15; %0314, 0319
I.time = datetime([cfg.ymd,0,0,dasc.seconds_since_midnight]);
I.mlon = dasc.maglon_110(id00:id01,1)+360;
I.mlat = dasc.maglat_110(1,id10:id11)';
I.glon = dasc.footlon_110(id00:id01,id10:id11);
I.glat = dasc.footlat_110(id00:id01,id10:id11);
I.flux = smoothdata2(dasc.Q(id00:id01,id10:id11),"Gaussian",sp);
I.energy = smoothdata2(dasc.E0(id00:id01,id10:id11),"Gaussian",sp);
I.pedersen = smoothdata2(dasc.SigP(id00:id01,id10:id11),"Gaussian",sp);
I.hall = smoothdata2(dasc.SigH(id00:id01,id10:id11),"Gaussian",sp);
% I.flux = dasc.Q(id00:id01,id10:id11);
% I.energy = dasc.E0(id00:id01,id10:id11);
% I.pedersen = dasc.SigP(id00:id01,id10:id11);
% I.hall = dasc.SigH(id00:id01,id10:id11);
I.red = dasc.red_rayleighs(id00:id01,id10:id11);
I.green = dasc.green_rayleighs(id00:id01,id10:id11);
I.blue = dasc.blue_rayleighs(id00:id01,id10:id11);

[MLON,MLAT] = ndgrid(I.mlon,I.mlat);
close all
figure('Position',[10 50 1500 800])
tiledlayout(2,2)
nexttile
pcolor(MLON,MLAT,I.flux); shading flat; colorbar
nexttile
pcolor(MLON,MLAT,I.energy); shading flat; colorbar
nexttile
pcolor(MLON,MLAT,I.pedersen); shading flat; colorbar
nexttile
pcolor(MLON,MLAT,I.hall); shading flat; colorbar

any_nans = [any(isnan(I.flux),'all')
any(isnan(I.energy),'all')
any(isnan(I.pedersen),'all')
any(isnan(I.hall),'all')];
if any(any_nans)
    error('Some NaNs left in image data')
end

%%
file = fullfile(direc,'ext','precipitation.h5');

h5make(file,'/Time/Year',year(I.time),'Year',type='int16')
h5make(file,'/Time/DOY',day(I.time,'dayofyear'),'Day of year',type='int16')
h5make(file,'/Time/Seconds',second(I.time,'secondofday'),'Seconds since midnight')
h5make(file,'/Time/Unix',posixtime(I.time),'Unix time')

h5make(file,'/Coordinates/Magnetic/Longitude',I.mlon,'Magnetic longitude' ...
    ,units='Degrees east (0, 360)',size='1 x Nmlon')
h5make(file,'/Coordinates/Magnetic/Latitude',I.mlat,'Magnetic latitude' ...
    ,units='Degrees north (-90, 90)',size='1 x Nmlat')
h5make(file,'/Coordinates/Geodetic/Longitude',I.glon,'Geodetic longitude' ...
    ,units='Degrees east (0, 360)',size='Nmlon x Nmlat',foot_alt='110 km')
h5make(file,'/Coordinates/Geodetic/Latitude',I.glat,'Geodetic latitude' ...
    ,units='Degrees north (-90, 90)',size='Nmlon x Nmlat',foot_alt='110 km')

h5make(file,'/Optical/Red',I.red,'Red filter',units='Rayleighs',size='Nmlon x Nmlat')
h5make(file,'/Optical/Green',I.green,'Green filter',units='Rayleighs',size='Nmlon x Nmlat')
h5make(file,'/Optical/Blue',I.blue,'Blue filter',units='Rayleighs',size='Nmlon x Nmlat')

h5make(file,'/Derived/Energy/Flux',I.flux,'Total precipitating energy flux' ...
    ,units='Milliwatts/meter^2',size='Nmlon x Nmlat')
h5make(file,'/Derived/Energy/Characteristic',I.energy,'Precipitating characteristic energy' ...
    ,units='Electronvolts',size='Nmlon x Nmlat')
h5make(file,'/Derived/Conductance/Pedersen',I.pedersen,'Pedersen conductance' ...
    ,units='Siemens',size='Nmlon x Nmlat')
h5make(file,'/Derived/Conductance/Hall',I.hall,'Hall conductance' ...
    ,units='Siemens',size='Nmlon x Nmlat')

%% swarm track data
file = fullfile(direc,'ext','tracks.h5');

for sat = 'ABC'
    try
        swarm_fn.(sat) = dir(fullfile(swarm_path,['*EFI',sat,'*2023',md,'*.h5'])).name;
    catch
        warning('No data found for Swarm %s.',sat)
    end
end

close all
figure
hold on
pcolor(MLON,MLAT,I.hall); shading flat

sats = cell2mat(fieldnames(swarm_fn))';
for sat = sats
    h5f = fullfile(swarm_path,swarm_fn.(sat));
    times = datetime(h5read(h5f,'/Timestamp'),'ConvertFrom','posixtime');
    mlon = h5read(h5f,'/MagneticLongitude');
    mlat = h5read(h5f,'/MagneticLatitude');
    glon = h5read(h5f,'/Longitude');
    glat = h5read(h5f,'/GeodeticLatitude');
    glon_foot = h5read(h5f,'/Longitude110km');
    glat_foot = h5read(h5f,'/GeodeticLatitude110km');
    mflow_east = h5read(h5f,'/ViMagE');
    mflow_north = h5read(h5f,'/ViMagN');
    gflow_east = h5read(h5f,'/ViE');
    gflow_north = h5read(h5f,'/ViN');
    fac = zeros(size(times)); % TBD
    
    buf = 0.5; % degrees mlat
    mflow_mag = sqrt(mflow_east.^2 + mflow_north.^2);
    max_flow = 6e3;
    max_diff_flow = 300;
    ok1 = mlat > (gem_mlat(1)-buf) & mlat < (gem_mlat(end) + buf);
    ok2 = mflow_mag < max_flow;
    ok3 = [abs(diff(mflow_mag)); nan] < max_diff_flow;
    ok = ok1 & ok2 & ok3;

    if any([sum(diff(ok)==1),sum(diff(ok)==-1)]>2)
        error('"ok" not consecutive.')
    end
    
    h5make(file,['/',sat,'/Time/Year'],year(times(ok)),'Year',type='int16')
    h5make(file,['/',sat,'/Time/DOY'],day(times(ok),'dayofyear'),'Day of year',type='int16')
    h5make(file,['/',sat,'/Time/Seconds'],second(times(ok),'secondofday'),'Seconds since midnight')
    h5make(file,['/',sat,'/Time/Unix'],posixtime(times(ok)),'Unix time')

    h5make(file,['/',sat,'/Coordinates/Magnetic/Longitude'],mlon(ok),'Magnetic longitude' ...
        ,units='Degrees east (0, 360)')
    h5make(file,['/',sat,'/Coordinates/Magnetic/Latitude'],mlat(ok),'Magnetic latitude' ...
        ,units='Degrees north (-90, 90)')
    h5make(file,['/',sat,'/Coordinates/Geodetic/Longitude'],glon(ok),'Geodetic longitude' ...
        ,units='Degrees east (0, 360)')
    h5make(file,['/',sat,'/Coordinates/Geodetic/Latitude'],glat(ok),'Geodetic latitude' ...
        ,units='Degrees north (-90, 90)')
    h5make(file,['/',sat,'/Coordinates/Geodetic/FootLongitude'],glon_foot(ok),'Footpointed Geodetic longitude' ...
        ,units='Degrees east (0, 360)',foot_alt='110 km')
    h5make(file,['/',sat,'/Coordinates/Geodetic/FootLatitude'],glat_foot(ok),'Footpointed Geodetic latitude' ...
        ,units='Degrees north (-90, 90)',foot_alt='110 km')

    h5make(file,['/',sat,'/Flow/Magnetic/East'],mflow_east(ok),'Magnetic eastward plasma flow',units='Meters/second')
    h5make(file,['/',sat,'/Flow/Magnetic/North'],mflow_north(ok),'Magnetic northward plasma flow',units='Meters/second')
    h5make(file,['/',sat,'/Flow/Geodetic/East'],gflow_east(ok),'Geodetic eastward plasma flow',units='Meters/second')
    h5make(file,['/',sat,'/Flow/Geodetic/North'],gflow_north(ok),'Geodetic northward plasma flow',units='Meters/second')
    h5make(file,['/',sat,'/Current/FieldAligned'],fac(ok),'Field aligned current',units='Microamperes/meter^2')

    scl = 1e-3;
    quiver(mlon(ok),mlat(ok),mflow_east(ok)*scl,mflow_north(ok)*scl,0,'.-r')
end

%% OLD
% buf = 1; % mlat
% ok = swarm.Mag_Lat > (I.mlat(1)-buf) & swarm.Mag_Lat < (I.mlat(end)+buf);
% ok = 386:686;

% UTtime = pad(string(num2str(swarm.UTtime')),10,'left','0');
% dt = datetime('today') - datetime(swarm_fn(20:27),'InputFormat','uuuuMMdd');
% A.times = datetime(UTtime,'InputFormat','HHmmss.SSS')' - dt;
% % A.times = datetime(swarm.UnixTime(ok),'ConvertFrom','posixtime');
% A.mlon = swarm.Mag_Lon(ok) + 360;
% A.mlat = swarm.Mag_Lat(ok);
% A.glon = swarm.Geo_Lon(ok) + 360;
% A.glat = swarm.Geo_Lat(ok);
% A.glon_foot = swarm.Footpoint_Lon(ok) + 360;
% A.glat_foot = swarm.Footpoint_Lat(ok);
% A.mflow_east = swarm.Flow_MagEast(ok);
% A.mflow_north = swarm.Flow_GeoNorth(ok);
% A.gflow_east = swarm.Flow_GeoEast(ok);
% A.gflow_north = swarm.Flow_MagNorth(ok);
% A.fac = swarm.FAC_micro(ok);

%%

% h5fn = 'plots\swop\SW_EXPT_EFIB_TCT02_20230304T100800.h5';
% mlon = h5read(h5fn,'/mlon');
% mlat = h5read(h5fn,'/mlat');
% v_geo_e = h5read(h5fn,'/v_geo_e');
% v_geo_n = h5read(h5fn,'/v_geo_n');
% v_mag_e = h5read(h5fn,'/v_mag_e');
% v_mag_n = h5read(h5fn,'/v_mag_n');

%% pfisr track data
% file_pfisr = 'data\20230319.001_lp_3min-fitcal-vvels_lat.h5';
% 
% B.mlon = h5read(file_pfisr,'/VvelsMagCoords/Longitude') + 360;
% B.mlat = h5read(file_pfisr,'/VvelsMagCoords/Latitude');
% mflow = h5read(file_pfisr,'/VvelsMagCoords/Velocity');
% mflow_err = h5read(file_pfisr,'/VvelsMagCoords/errVelocity');
% gflow = h5read(file_pfisr,'/VvelsGeoCoords/Velocity');
% alt = h5read(file_pfisr,'/VvelsGeoCoords/Altitude');
% 
% tid = 34;
% ok1 = mean(mflow_err(:,:,tid))' < 300;
% ok2 = B.mlat > (I.mlat(1)-buf) & B.mlat < (I.mlat(end)+buf);
% ok = ok1 & ok2;
% l_vvels = length(B.mlon(ok));
% 
% B.times = h5read(file_pfisr,'/Time/UnixTime');
% 
% B.time = datetime(B.times(:,tid)','ConvertFrom','posixtime');
% B.mlon = B.mlon(ok) + 1e-5*(0:l_vvels-1)';
% B.mlat = B.mlat(ok);
% B.mflow_east = squeeze(mflow(1,ok,tid))'; % see Laundal & Richmond (2016), Eq. 57-59
% B.mflow_north = -squeeze(mflow(2,ok,tid))';

%% save track file
% file = fullfile(direc,'ext','tracks_pfisr-only3.h5');
% 
% h5make(file,'/NumTracks',1,'Number of tracks',type='int16')

% track A
% h5make(file,'/A/Time/Year',year(A.times),'Year',type='int16')
% h5make(file,'/A/Time/DOY',day(A.times,'dayofyear'),'Day of year',type='int16')
% h5make(file,'/A/Time/Seconds',second(A.times,'secondofday'),'Seconds since midnight')
% h5make(file,'/A/Time/Unix',posixtime(A.times),'Unix time')
% 
% h5make(file,'/A/Coordinates/Magnetic/Longitude',A.mlon,'Magnetic longitude' ...
%     ,units='Degrees east (0, 360)')
% h5make(file,'/A/Coordinates/Magnetic/Latitude',A.mlat,'Magnetic latitude' ...
%     ,units='Degrees north (-90, 90)')
% h5make(file,'/A/Coordinates/Geodetic/Longitude',A.glon,'Geodetic longitude' ...
%     ,units='Degrees east (0, 360)')
% h5make(file,'/A/Coordinates/Geodetic/Latitude',A.glat,'Geodetic latitude' ...
%     ,units='Degrees north (-90, 90)')
% h5make(file,'/A/Coordinates/Geodetic/FootLongitude',A.glon_foot,'Footpointed Geodetic longitude' ...
%     ,units='Degrees east (0, 360)',foot_alt='110 km')
% h5make(file,'/A/Coordinates/Geodetic/FootLatitude',A.glat_foot,'Footpointed Geodetic latitude' ...
%     ,units='Degrees north (-90, 90)',foot_alt='110 km')
% 
% h5make(file,'/A/Flow/Magnetic/East',A.mflow_east,'Magnetic eastward plasma flow',units='Meters/second')
% h5make(file,'/A/Flow/Magnetic/North',A.mflow_north,'Magnetic northward plasma flow',units='Meters/second')
% h5make(file,'/A/Flow/Geodetic/East',A.gflow_east,'Geodetic eastward plasma flow',units='Meters/second')
% h5make(file,'/A/Flow/Geodetic/North',A.gflow_north,'Geodetic northward plasma flow',units='Meters/second')
% h5make(file,'/A/Current/FieldAligned',A.fac,'Field aligned current',units='Microamperes/meter^2')

% track B
% h5make(file,'/B/Time/Year',year(B.time),'Year' ...
%     ,type='int16',size='1 x 2 (start, end)')
% h5make(file,'/B/Time/DOY',day(B.time,'dayofyear'),'Day of year' ...
%     ,type='int16',size='1 x 2 (start, end)')
% h5make(file,'/B/Time/Seconds',second(B.time,'secondofday'),'Seconds since midnight' ...
%     ,size='1 x 2 (start, end)')
% h5make(file,'/B/Time/Unix',posixtime(B.time),'Unix time' ...
%     ,size='1 x 2 (start, end)')
% 
% h5make(file,'/B/Coordinates/Magnetic/Longitude',B.mlon,'Magnetic longitude' ...
%     ,units='Degrees east (0, 360)')
% h5make(file,'/B/Coordinates/Magnetic/Latitude',B.mlat,'Magnetic latitude' ...
%     ,units='Degrees north (-90, 90)')
% 
% h5make(file,'/B/Flow/Magnetic/East',B.mflow_east,'Magnetic eastward plasma flow' ...
%     ,units='Meters/second')
% h5make(file,'/B/Flow/Magnetic/North',B.mflow_north,'Magnetic northward plasma flow' ...
%     ,units='Meters/second')
