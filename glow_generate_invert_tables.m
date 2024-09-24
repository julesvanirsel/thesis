glow_direc = fullfile('data', 'paper2', 'dasc_data', 'glow');
event_fn = fullfile('data', 'paper2', 'event_data.txt');
site_lat = 65.12;
site_lon = 212.57;
lines = readlines(event_fn)';
for l = lines(2:end-1)
    data = strsplit(l);
    date = datetime(data(2), 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss''Z');
    jules.glow.generate_invert_table(date, site_lat, site_lon, out_direc=glow_direc)
end