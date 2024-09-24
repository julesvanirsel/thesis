import numpy as np
import matplotlib.pyplot as plt
from asispectralinversion import preprocessing, inversion
import apexpy
from datetime import datetime, timedelta, timezone
import os
import requests
import imageio.v3 as iio
import h5py


def main(direc, date, grid_size_lat, grid_size_lon, num_coadd=3, time_window=15, num_shifts=30, background_method='patches', debug=False, dec = 2):

    # init
    if date.tzinfo == None:
        raise Exception('Date requires timezone. E.g. date = datetime(..., tzinfo=timezone.utc).')
    direc = make_valid_path(direc)
    A = apexpy.Apex(date=date)

    # load skymap
    print_progress('LOADING SKYMAP')
    skymap_red, skymap_grn, skymap_blu, site_lat, site_lon, site_alt = load_skymap(direc)
    site_mag_lat, _ = A.convert(site_lat, site_lon, 'geo', 'apex', height=site_alt)
    print('DONE')

    # load glow inversion lookup tables
    print_progress('LOADING INVERSION LOOKUP TABLES')
    glow_path = os.path.join(direc, 'glow', date.strftime('%y%j_') + str(sod(date))) + os.sep
    if debug:
        print('\nGLOW PATH = {}'.format(glow_path))
    invert_table = inversion.load_lookup_tables_directory(glow_path, site_mag_lat)
    print('DONE')

    # download imagery
    img_times = download_imagery(direc, date, time_window)

    # coadd images
    print_progress('COADDING IMAGERY')
    coadd_times_red, time_red = get_coadd_times(date, img_times['630'], num_coadd)
    coadd_times_grn, time_grn = get_coadd_times(date, img_times['558'], num_coadd)
    coadd_times_blu, time_blu = get_coadd_times(date, img_times['428'], num_coadd)
    img_red = load_image(direc, coadd_times_red, 630)
    img_grn = load_image(direc, coadd_times_grn, 558)
    img_blu = load_image(direc, coadd_times_blu, 428)
    print('DONE')

    # choose inversion grid
    lims_lat, lims_lon = A.convert(
        [np.nanmin(skymap_blu[0]), np.nanmax(skymap_blu[0])], [np.nanmin(skymap_blu[1]), np.nanmax(skymap_blu[1])], 'geo', 'apex', height=110)
    interp_lat = np.linspace(lims_lat[0], lims_lat[1], grid_size_lat*dec)
    interp_lon = np.linspace(lims_lon[0], lims_lon[1], grid_size_lon*dec)

    # wavelet denoise the images
    print_progress('WAVELET DENOISING')
    img_denoise_red, _, geod_lon, geod_lat, _, _, bg_red, sigma_red = preprocessing.wavelet_denoise_resample(
        img_red, date, skymap_red[1], skymap_red[0], interp_lon, interp_lat, 110, nshifts=num_shifts, background_method=background_method)
    img_denoise_grn, _, _, _, _, _, bg_grn, sigma_grn = preprocessing.wavelet_denoise_resample(
        img_grn, date, skymap_grn[1], skymap_grn[0], interp_lon, interp_lat, 110, nshifts=num_shifts, background_method=background_method)
    img_denoise_blu, _, _, _, _, _, bg_blu, sigma_blu = preprocessing.wavelet_denoise_resample(
        img_blu, date, skymap_blu[1], skymap_blu[0], interp_lon, interp_lat, 110, nshifts=num_shifts, background_method=background_method)
    print('DONE')

    # convert to raylaighs
    ray_red, ray_grn, ray_blu = preprocessing.to_rayleighs(img_denoise_red, img_denoise_grn, img_denoise_blu, bg_red, bg_grn, bg_blu)

    # remove non-positive values
    ray_red[np.where(ray_red <= 0)] = np.nan
    ray_grn[np.where(ray_grn <= 0)] = np.nan
    ray_blu[np.where(ray_blu <= 0)] = np.nan

    # remove values where one or more color is missing
    nan_ids = np.where(np.isnan(ray_red + ray_grn + ray_blu))
    ray_red[nan_ids] = np.nan
    ray_grn[nan_ids] = np.nan
    ray_blu[nan_ids] = np.nan

    # decimate down to final size
    geod_lat= geod_lat[::dec, ::dec]
    geod_lon= geod_lon[::dec, ::dec]
    ray_red_dec = ray_red[::dec, ::dec]
    ray_grn_dec = ray_grn[::dec, ::dec]
    ray_blu_dec = ray_blu[::dec, ::dec]
    img_denoise_red = img_denoise_red[::dec, ::dec]
    img_denoise_grn = img_denoise_grn[::dec, ::dec]
    img_denoise_blu = img_denoise_blu[::dec, ::dec]

    # invert imagery
    print_progress('INVERTING IMAGERY')
    Q, E0, Q_bound_lo, Q_bound_hi, E0_bound_lo, E0_bound_hi = inversion.calculate_E0_Q_v2(ray_red_dec, ray_grn_dec, ray_blu_dec, invert_table, minE0=150, generous=False)
    SigP, SigH = inversion.calculate_Sig(Q, E0, invert_table, generous=False)
    print('DONE')

    # saving inversions
    print_progress('SAVING INVERSION DATA')
    optical_times = [time_red.timestamp(), time_grn.timestamp(), time_blu.timestamp()]
    optical_coadd_times = [[t.timestamp() for t in coadd_times_red], [t.timestamp() for t in coadd_times_grn], [t.timestamp() for t in coadd_times_blu]]
    h5f_out = h5py.File(os.path.join(direc, 'INV_PKR_' + date.strftime('%Y%m%d_') + str(sod(date)) + '.h5'), 'w')

    h5f_out['/'].attrs['GLOBAL DATUM'] = 'WGS84'
    h5f_out['/'].attrs['APEX VERSION'] = apexpy.__version__
    h5f_out['/'].attrs['GLOW PATH'] = glow_path
    h5f_out['/'].attrs['GRID SIZE X'] = grid_size_lon
    h5f_out['/'].attrs['GRID SIZE Y'] = grid_size_lat
    h5f_out['/'].attrs['IMAGE DOWNLOAD TIME RANGE'] = str(2*time_window) + ' seconds'
    h5f_out['/'].attrs['BACKGROUND METHOD'] = background_method
    h5f_out['/'].attrs['Nx'] = 'Number of cells approximately eastward'
    h5f_out['/'].attrs['Ny'] = 'Number of cells approximately nothward'
    h5f_out['/'].attrs['Nxi'] = 'Number of image pixels approximately nothward'
    h5f_out['/'].attrs['Nyi'] = 'Number of image pixels approximately nothward'
    h5f_out['/'].attrs['Nc'] = 'Number of coadds'

    h5_write(h5f_out, 'Coordinates/Latitude', geod_lat, 'float64', 'degrees north', 'Geodetic latitude', 'Nx x Ny')
    h5_write(h5f_out, 'Coordinates/Longitude', np.mod(geod_lon, 360), 'float64', 'degrees east', 'Geodetic longitude', 'Nx x Ny')
    h5_write(h5f_out, 'Coordinates/Altitude', 110e3, 'float64', 'meters', 'Geodetic altitude', 'Scalar')

    h5_write(h5f_out, 'Time/Unix', np.mean(optical_times), 'float64', 'seconds', 'Unix time of inversion', 'Scalar')
    h5_write(h5f_out, 'Time/Year', date.year, 'int16', 'counts', 'Year', 'Scalar')
    h5_write(h5f_out, 'Time/DOY', doy(date), 'int16', 'counts', 'Day of year', 'Scalar')
    h5_write(h5f_out, 'Time/SOD', sod(date), 'float64', 'seconds', 'Second of day', 'Scalar')

    h5_write(h5f_out, 'Derived/Energy/Flux', Q, 'float64', 'milliwatts/meter^2', 'Precipitating electron energy flux', 'Nx x Ny')
    h5_write(h5f_out, 'Derived/Energy/FluxBounds', [Q_bound_lo, Q_bound_hi], 'float64', 'milliwatts/meter^2', 'Precipitating electron energy flux bounds', '2 (low, high) x Nx x Ny')
    h5_write(h5f_out, 'Derived/Energy/Characteristic', E0, 'float64', 'electronvolts', 'Precipitating electron characteristic energy', 'Nx x Ny')
    h5_write(h5f_out, 'Derived/Energy/CharacteristicBounds', [E0_bound_lo, E0_bound_hi], 'float64', 'electronvolts', 'Precipitating electron characteristic energy bounds', '2 (low, high) x Nx x Ny')
    h5_write(h5f_out, 'Derived/Conductance/Pedersen', SigP, 'float64', 'siemens', 'Pedersen Conductance', 'Nx x Ny')
    h5_write(h5f_out, 'Derived/Conductance/Hall', SigH, 'float64', 'siemens', 'Hall Conductance', 'Nx x Ny')

    h5_write(h5f_out, 'Optical/Skymaps/Latitude', [skymap_red[0], skymap_grn[0], skymap_blu[0]], 'float64', 'degrees north', 'Geodetic latitude', '3 (red, green, blue) x Nxi x Nyi')
    h5_write(h5f_out, 'Optical/Skymaps/Longitude', np.mod([skymap_red[1], skymap_grn[1], skymap_blu[1]], 360), 'float64', 'degrees east', 'Geodetic longitude', '3 (red, green, blue) x Nxi x Nyi')
    h5_write(h5f_out, 'Optical/Skymaps/Altitude', 110e3, 'float64', 'meters', 'Geodetic altitude', 'Scalar')
    h5_write(h5f_out, 'Optical/Red', ray_red_dec, 'float64', 'rayleighs', 'Red digital all-sky image', 'Nx x Ny')
    h5_write(h5f_out, 'Optical/Green', ray_grn_dec, 'float64', 'rayleighs', 'Green digital all-sky image', 'Nx x Ny')
    h5_write(h5f_out, 'Optical/Blue', ray_blu_dec, 'float64', 'rayleighs', 'Blue digital all-sky image', 'Nx x Ny')
    h5_write(h5f_out, 'Optical/Times', optical_times, 'float64', 'seconds', 'Unix times of imagery', '3 (red, green, blue)')
    h5_write(h5f_out, 'Optical/Wavelengths', [630e-9, 558e-9, 428e-9], 'float64', 'meters', 'Photon wavelengths of images', '3 (red, green, blue)')
    h5_write(h5f_out, 'Optical/Preprocessing/NumShifts', num_shifts, 'int16', 'counts', '"max_shifts" from skimage.restoration.cycle_spin', 'Scalar')
    h5_write(h5f_out, 'Optical/Preprocessing/NumCoadd', num_coadd, 'int16', 'counts', 'Number of coadded images per color', 'Scalar')
    h5_write(h5f_out, 'Optical/Preprocessing/Coadded', [img_red, img_grn, img_blu], 'float64', 'counts', 'Coadded images', '3 (red, green, blue) x Nxi x Nyi')
    h5_write(h5f_out, 'Optical/Preprocessing/CoaddTimes', optical_coadd_times, 'float64', 'seconds', 'Unix times of coadded images', '3 (red, green, blue) x Nc')
    h5_write(h5f_out, 'Optical/Preprocessing/Denoised', [img_denoise_red, img_denoise_grn, img_denoise_blu], 'float64', 'counts', 'Denoised images', '3 (red, green, blue) x Nx x Ny')
    h5_write(h5f_out, 'Optical/Preprocessing/Sigma', [sigma_red, sigma_grn, sigma_blu], 'float64', 'counts', 'Gaussian standard deviations of images backgrounds', '3 (red, green, blue)')
    h5_write(h5f_out, 'Optical/Preprocessing/Background', [bg_red, bg_grn, bg_blu], 'float64', 'counts', 'Background values used in denoising', '3 (red, green, blue)')

    h5_write(h5f_out, 'Site/Latitude', site_lat, 'float64', 'degrees north', 'Geodetic latitude of site', 'Scalar')
    h5_write(h5f_out, 'Site/Longitude', np.mod(site_lon, 360), 'float64', 'degrees east', 'Geodetic longitude of site', 'Scalar')
    h5_write(h5f_out, 'Site/Altitude', site_alt, 'float64', 'meters', 'Geodetic altitude of site', 'Scalar')
    h5_write(h5f_out, 'Site/MagneticLatitude', site_mag_lat, 'float64', 'degrees north', 'Apex magnetic latitude of site', 'Scalar')

    print('DONE')


def get_coadd_times(date, times, num):
    if len(times) == num:
        coadd_times = sorted(times)
        coadd_time = coadd_times[0] + timedelta(0, (coadd_times[-1]-coadd_times[0]).total_seconds() / 2)
        return coadd_times, coadd_time
    elif len(times) < num:
        raise Exception('Not enough frames to coadd.')
    j = 0
    coadd_times = [-1]*num
    deltas = [abs((t-date).total_seconds()) for t in times]
    for d in sorted(deltas)[:num]:
        for i in range(len(deltas)):
            if d == deltas[i]:
                coadd_times[j] = times[i]
                j += 1
                deltas[i] = -1
                break
    coadd_times = sorted(coadd_times)
    coadd_time = coadd_times[0] + timedelta(0, (coadd_times[-1]-coadd_times[0]).total_seconds() / 2)
    return coadd_times, coadd_time


def download_imagery(direc, date, window, debug=False):
    dasc_root = 'http://optics.gi.alaska.edu/amisr_archive/PKR/DASC/PNG/'
    times = {'428':[], '558':[], '630':[]}
    tot = (window*2 + 1)*3
    pct = 0
    for s in range(-window, window + 1):
        for c in [428, 558, 630]:
            t = date + timedelta(0, s)
            url = dasc_root + t.strftime('%Y/%Y%m%d/%H/')
            filename = t.strftime('PKR_%Y%m%d_%H%M%S_{}.png'.format(str(c).rjust(4, '0')))
            out_path = make_valid_path(os.path.join(direc, 'imagery', t.strftime('%Y%m%d')))
            out_path = os.path.join(out_path, filename)
            if os.path.exists(out_path):
                pct += 1
                if not(debug):
                    load_bar('DOWNLOADING IMAGERY', pct, tot)
                times[str(c)].append(t)
                continue
            image = requests.get(url + filename)
            if image.status_code == 200:
                times[str(c)].append(t)
                with open(out_path, 'wb') as handler:
                    handler.write(image.content)
            elif (image.status_code == 404) & debug:
                print('\n' + url + filename + ' NOT FOUND.')
            elif debug:
                print('\nREQUEST RETURNED STATUS CODE ' + str(image.status_code) + '.')
            pct += 1
            if not(debug):
                load_bar('DOWNLOADING IMAGERY', pct, tot)
    print('DONE')
    return times


def load_image(direc, times, color):
    n = 1
    for t in times:
        filename = t.strftime('PKR_%Y%m%d_%H%M%S_{}.png'.format(str(color).rjust(4, '0')))
        out_path = os.path.join(direc, 'imagery', t.strftime('%Y%m%d'), filename)
        img_now = iio.imread(out_path)
        if n > 1:
            img_avg = (img_now + img_prv*(n-1))/n
        else:
            img_avg = img_now
        img_prv = img_avg
        n += 1
    return img_avg


def load_skymap(direc, debug=False):
    h5f = h5py.File(os.path.join(direc, 'imagery', 'skymap.h5'), 'r')
    site_lat = h5f['Site/Latitude'][0]
    site_lon = h5f['Site/Longitude'][0]
    site_alt = h5f['Site/Altitude'][0] / 1e3
    lat = h5f['Footpointing/Apex/Latitude']
    lon = h5f['Footpointing/Apex/Longitude']
    alt = h5f['Footpointing/Apex/Altitude']
    skymap_red = [lat[2, :, 0:-1], lon[2, :, 0:-1]] # map deliberately assymetric to id axes
    skymap_grn = [lat[1, :, 0:-1], lon[1, :, 0:-1]]
    skymap_blu = [lat[0, :, 0:-1], lon[0, :, 0:-1]]
    if debug:
        print('SKYMAP SOURCE ALTITUDES: {} METERS'.format(np.transpose(alt)))
    return skymap_red, skymap_grn, skymap_blu, site_lat, site_lon, site_alt


def geod_apex_angle(geod_lat, geod_lon, geod_alt, date):
    A = apexpy.Apex(date=date)
    _, _, _, _, _, _, _, _, _, _, e2, _ = A.basevectors_apex(geod_lat, geod_lon, geod_alt) # geodetic east, north, up
    if geod_lat > 0:
        emN = -e2 # e2 points equatorward
    else:
        emN = e2
    angle = np.arccos(emN[1])
    print('Apex magnetic east, north, up vector is ({:.4f}, {:.4f}, {:.4f}).'.format(emN[0], emN[1], emN[2]))
    print('The angle between magnetic and geodetic north is {:.4f}° east of north.'.format(angle*180/np.pi))
    return angle


def make_valid_path(path):
    path = os.path.abspath(path)
    if not(os.path.exists(path)):
        try:
            os.mkdir(path)
        except FileNotFoundError:
            try:
                os.mkdir(os.path.split(path)[0])
                os.mkdir(path)
            except FileNotFoundError:
                print('COULD NOT CREATE PATH: ' + path)
    return path


def h5_write(file, name, data, type, units, description, size):
    ds = file.create_dataset(name, data=np.squeeze(np.array(data)), dtype=type)
    ds.attrs['DESCRIPTION'] = description
    ds.attrs['UNITS'] = units
    ds.attrs['SIZE'] = size


def load_bar(msg, pct, tot, wdth=75):
    print((msg + ': ' + '.' * round((wdth - len(msg) - 2) * pct / tot)).ljust(wdth), end=' ', flush=True)
    print('{} / {} ({:.0f}%)'.format(pct, tot, 100 * pct / tot) + '\r', end='')


def print_progress(name, wdth=75):
    print((name + ' ').ljust(wdth, '.'), end=' ', flush=True)  


def doy(date):
    return (date - datetime(date.year, 1, 1, tzinfo=date.tzinfo)).days + 1


def sod(date):
    return (date - datetime(date.year, 1, 1, tzinfo=date.tzinfo)).seconds


if __name__ == '__main__':
    direc = os.path.join('data', 'paper2', 'dasc_data')
    grid_size_lat = 512+128
    grid_size_lon = 512
    with open(os.path.join('data', 'paper2', 'event_data.txt'), 'r') as f:
        lines = f.readlines()
        for l in lines[1:2]:
            date_str = l.split()[1].replace('Z','+0000')
            date = datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S%z')
            print(date)
            main(direc, date, grid_size_lat, grid_size_lon)