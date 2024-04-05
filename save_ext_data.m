direc = '\\dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\test_run';
h5make = @jules.tools.h5make;

% imagery data
load('data\SW_EXPT_EFIB_TCT02_20230319T081600_sky_inversion.mat')
dasc = event;
clear('event')
ids0 = [54,88]; % mlat
ids1 = [27,94]; % mlon

I.time = datetime(2023,3,19,0,0,30164+60);
I.mlon = dasc.maglon_110(1,ids1(1):ids1(2))+360;
I.mlat = dasc.maglat_110(ids0(1):ids0(2),1)';
I.glon = dasc.footlon_110(ids0(1):ids0(2),ids1(1):ids1(2))';
I.glat = dasc.footlat_110(ids0(1):ids0(2),ids1(1):ids1(2))';
I.flux = dasc.Q(ids0(1):ids0(2),ids1(1):ids1(2))';
I.energy = dasc.E0(ids0(1):ids0(2),ids1(1):ids1(2))';
I.pedersen = dasc.SigP(ids0(1):ids0(2),ids1(1):ids1(2))';
I.hall = nan(size(I.pedersen));
I.red = dasc.red_rayleighs(ids0(1):ids0(2),ids1(1):ids1(2))';
I.green = dasc.green_rayleighs(ids0(1):ids0(2),ids1(1):ids1(2))';
I.blue = dasc.blue_rayleighs(ids0(1):ids0(2),ids1(1):ids1(2))';

%% save precipitation file
file = fullfile(direc,'ext','precipitation.h5');

h5make(file,'/Time/Year',year(I.time),'Year',type='int16')
h5make(file,'/Time/DOY',day(I.time,'dayofyear'),'Day of year',type='int16')
h5make(file,'/Time/Seconds',second(I.time,'secondofday'),'Seconds since midnight')
h5make(file,'/Time/Unix',posixtime(I.time),'Unix time')

h5make(file,'/Coordinates/Magnetic/Longitude',I.mlon,'Magnetic longitude' ...
    ,units='Degrees east (0, 360)',size='1 x Nmlon')
h5make(file,'/Coordinates/Magnetic/Latitude',I.mlat,'Magnetic latitude' ...
    ,units='Degrees north (-90, 90)',size='1 x Nmlat')
h5make(file,'/Coordinates/Geographic/Longitude',I.glon,'Geographic longitude' ...
    ,units='Degrees east (0, 360)',size='Nmlon x Nmlat',foot_alt='110 km')
h5make(file,'/Coordinates/Geographic/Latitude',I.glat,'Geographic latitude' ...
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
load('data\SW_EXPT_EFIB_TCT02_20230319T081600_alongtrackincluded_used.mat')
swarm = SW_EXPT_EFIB_TCT02_20230319T081600;
clear('SW_EXPT_EFIB_TCT02_20230319T081600')

buf = 1; % mlat
ok = swarm.Mag_Lat > (I.mlat(1)-buf) & swarm.Mag_Lat < (I.mlat(end)+buf);

A.times = datetime(swarm.UnixTime(ok),'ConvertFrom','posixtime');
A.mlon = swarm.Mag_Lon(ok) + 360;
A.mlat = swarm.Mag_Lat(ok);
A.glon = swarm.Geo_Lon(ok) + 360;
A.glat = swarm.Geo_Lat(ok);
A.glon_foot = swarm.Footpoint_Lon(ok) + 360;
A.glat_foot = swarm.Footpoint_Lat(ok);
A.mflow_east = swarm.Flow_MagEast(ok);
A.mflow_north = swarm.Flow_GeoNorth(ok);
A.gflow_east = swarm.Flow_GeoEast(ok);
A.gflow_north = swarm.Flow_MagNorth(ok);
A.fac = swarm.FAC_micro(ok);

%% pfisr track data
file_pfisr = 'data\20230319.001_lp_3min-fitcal-vvels_lat.h5';

B.mlon = h5read(file_pfisr,'/VvelsMagCoords/Longitude')+360;
B.mlat = h5read(file_pfisr,'/VvelsMagCoords/Latitude');
mflow = h5read(file_pfisr,'/VvelsMagCoords/Velocity');
mflow_err = h5read(file_pfisr,'/VvelsMagCoords/errVelocity');
gflow = h5read(file_pfisr,'/VvelsGeoCoords/Velocity');
alt = h5read(file_pfisr,'/VvelsGeoCoords/Altitude');

tid = 34;
ok1 = mean(mflow_err(:,:,tid))' < 300;
ok2 = B.mlat > (I.mlat(1)-buf) & B.mlat < (I.mlat(end)+buf);
ok = ok1 & ok2;
l_vvels = length(B.mlon(ok));

B.times = h5read(file_pfisr,'/Time/UnixTime');

B.time = datetime(B.times(:,tid)','ConvertFrom','posixtime');
B.mlon = B.mlon(ok) + 1e-5*(0:l_vvels-1)';
B.mlat = B.mlat(ok);
B.mflow_east = squeeze(mflow(1,ok,tid))'; % see Laundal & Richmond (2016), Eq. 57-59
B.mflow_north = -squeeze(mflow(2,ok,tid))';

%% save track file
file = fullfile(direc,'ext','tracks.h5');

h5make(file,'/NumTracks',2,'Number of tracks',type='int16')

% track A
h5make(file,'/A/Time/Year',year(A.times),'Year',type='int16')
h5make(file,'/A/Time/DOY',day(A.times,'dayofyear'),'Day of year',type='int16')
h5make(file,'/A/Time/Seconds',second(A.times,'secondofday'),'Seconds since midnight')
h5make(file,'/A/Time/Unix',posixtime(A.times),'Unix time')

h5make(file,'/A/Coordinates/Magnetic/Longitude',A.mlon,'Magnetic longitude' ...
    ,units='Degrees east (0, 360)')
h5make(file,'/A/Coordinates/Magnetic/Latitude',A.mlat,'Magnetic latitude' ...
    ,units='Degrees north (-90, 90)')
h5make(file,'/A/Coordinates/Geographic/Longitude',A.glon,'Geographic longitude' ...
    ,units='Degrees east (0, 360)')
h5make(file,'/A/Coordinates/Geographic/Latitude',A.glat,'Geographic latitude' ...
    ,units='Degrees north (-90, 90)')
h5make(file,'/A/Coordinates/Geographic/FootLongitude',A.glon_foot,'Footpointed geographic longitude' ...
    ,units='Degrees east (0, 360)',foot_alt='110 km')
h5make(file,'/A/Coordinates/Geographic/FootLatitude',A.glat_foot,'Footpointed geographic latitude' ...
    ,units='Degrees north (-90, 90)',foot_alt='110 km')

h5make(file,'/A/Flow/Magnetic/East',A.mflow_east,'Magnetic eastward plasma flow',units='Meters/second')
h5make(file,'/A/Flow/Magnetic/North',A.mflow_north,'Magnetic northward plasma flow',units='Meters/second')
h5make(file,'/A/Flow/Geographic/East',A.gflow_east,'Geographic eastward plasma flow',units='Meters/second')
h5make(file,'/A/Flow/Geographic/North',A.gflow_north,'Geographic northward plasma flow',units='Meters/second')
h5make(file,'/A/Current/FieldAligned',A.fac,'Field aligned current',units='Microamperes/meter^2')

% track B
h5make(file,'/B/Time/Year',year(B.time),'Year' ...
    ,type='int16',size='1 x 2 (start, end)')
h5make(file,'/B/Time/DOY',day(B.time,'dayofyear'),'Day of year' ...
    ,type='int16',size='1 x 2 (start, end)')
h5make(file,'/B/Time/Seconds',second(B.time,'secondofday'),'Seconds since midnight' ...
    ,size='1 x 2 (start, end)')
h5make(file,'/B/Time/Unix',posixtime(B.time),'Unix time' ...
    ,size='1 x 2 (start, end)')

h5make(file,'/B/Coordinates/Magnetic/Longitude',B.mlon,'Magnetic longitude' ...
    ,units='Degrees east (0, 360)')
h5make(file,'/B/Coordinates/Magnetic/Latitude',B.mlat,'Magnetic latitude' ...
    ,units='Degrees north (-90, 90)')

h5make(file,'/B/Flow/Magnetic/East',B.mflow_east,'Magnetic eastward plasma flow' ...
    ,units='Meters/second')
h5make(file,'/B/Flow/Magnetic/North',B.mflow_north,'Magnetic northward plasma flow' ...
    ,units='Meters/second')
