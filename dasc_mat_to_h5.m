matfn = 'data\paper2\dasc_data\imagery\skymap.mat';
load(matfn)
h5fn = [matfn(1:end-3), 'h5'];

% add row of nans to make assymetric grid
az = [az'; nan(1, 512)];
el = [el'; nan(1, 512)];

vlat107 = [vertical_footpointing.('107km').lat'; nan(1, 512)];
vlat110 = [vertical_footpointing.('110km').lat'; nan(1, 512)];
vlat180 = [vertical_footpointing.('180km').lat'; nan(1, 512)];
vlon107 = [vertical_footpointing.('107km').lon'; nan(1, 512)];
vlon110 = [vertical_footpointing.('110km').lon'; nan(1, 512)];
vlon180 = [vertical_footpointing.('180km').lon'; nan(1, 512)];
vlat = reshape([vlat107, vlat110, vlat180],[513,512,3]);
vlon = reshape([vlon107, vlon110, vlon180],[513,512,3]);
vlon = vlon + 360;

mlat107 = [magnetic_footpointing.('107km').lat'; nan(1, 512)];
mlat110 = [magnetic_footpointing.('110km').lat'; nan(1, 512)];
mlat180 = [magnetic_footpointing.('180km').lat'; nan(1, 512)];
mlon107 = [magnetic_footpointing.('107km').lon'; nan(1, 512)];
mlon110 = [magnetic_footpointing.('110km').lon'; nan(1, 512)];
mlon180 = [magnetic_footpointing.('180km').lon'; nan(1, 512)];
mlat = reshape([mlat107, mlat110, mlat180],[513,512,3]);
mlon = reshape([mlon107, mlon110, mlon180],[513,512,3]);
mlon = mlon + 360;

h5make(h5fn, 'Azimuth', az, 'double', 'degrees', 'Poker Flat DASC Azimuth', 'Ny x Nx')
h5make(h5fn, 'Elevation', el, 'double', 'degrees', 'Poker Flat DASC Elevation', 'Ny x Nx')
h5make(h5fn, 'Footpointing/Vertical/Latitude', vlat, 'double', 'degrees north', 'Linearly footpointed geodetic latitude', '3 x Ny x Nx (blue, green, red)')
h5make(h5fn, 'Footpointing/Vertical/Longitude', vlon, 'double', 'degrees east', 'Linearly footpointed longitude', '3 x Ny x Nx (blue, green, red)')
h5make(h5fn, 'Footpointing/Vertical/Altitude', [107, 110, 180]*1e3, 'double', 'meters', 'Geodetic altitude', '3 x 1 (blue, green, red)')
h5make(h5fn, 'Footpointing/Apex/Latitude', mlat, 'double', 'degrees north', 'Apex footpointed geodetic latitude', '3 x Ny x Nx (blue, green, red)')
h5make(h5fn, 'Footpointing/Apex/Longitude', mlon, 'double', 'degrees east', 'Apex footpointed longitude', '3 x Ny x Nx (blue, green, red)')
h5make(h5fn, 'Footpointing/Apex/Altitude', [107, 110, 180]*1e3, 'double', 'meters', 'Footpointing source geodetic altitude', '3 x 1 (blue, green, red)')
h5make(h5fn, 'Site/Latitude', 65.12, 'double', 'degrees north', 'Geodetic latitude of DASC site.', '1 x 1')
h5make(h5fn, 'Site/Longitude', 212.57, 'double', 'degrees east', 'Longitude of DASC site.', '1 x 1')
h5make(h5fn, 'Site/Altitude', 213, 'double', 'meters', 'Geodetic altitude of DASC site.', '1 x 1')

h5writeatt(h5fn, '/', 'Global datum', 'WGS84')
h5writeatt(h5fn, '/', 'Footpointing destination geodetic altitude', '110e3 meters')

function h5make(pth, var, dat, typ, unt, des, siz)
    % dat = permute(dat, length(size(dat)):-1:1);
    h5create(pth, ['/',var], size(dat), 'Datatype', typ)
    h5write(pth, ['/',var], dat)
    h5writeatt(pth, ['/',var], 'Units', unt)
    h5writeatt(pth, ['/',var], 'Description', des)
    h5writeatt(pth, ['/',var], 'Size', siz)
end

