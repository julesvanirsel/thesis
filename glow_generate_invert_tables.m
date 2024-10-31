glow_direc = fullfile('data', 'paper2', 'dasc_data', 'glow');
event_fn = fullfile('data', 'paper2', 'event_data.txt');
maph5 = fullfile('data', 'paper2', 'dasc_data', 'imagery', 'skymap.h5');
site_lat = h5read(maph5, '/Site/Latitude');
site_lon = h5read(maph5, '/Site/Longitude');
lines = readlines(event_fn)';

done_dates = NaT(1, length(lines)-2);
li0 = input('Start event ID: ');
do_acc = input('Accelerated Maxwellian? (0/1): ' );
for li = li0:length(lines)-2
    l = lines(li+1);
    data = strsplit(l);
    date = datetime(data(2), 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss''Z');
    te = str2double(data(13));
    fprintf('Source region thermal energy = %i eV\n', round(te))
    if any(done_dates == date)
        continue
    else
        done_dates(li) = date;
    end
    jules.glow.generate_invert_table(date, site_lat, site_lon, ...
        out_direc=glow_direc, do_acc=do_acc, thermal_energy_ev=te)
end