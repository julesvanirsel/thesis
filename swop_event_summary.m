direc = fullfile('data','paper2');
colorcet = @jules.tools.colorcet;
colors = [["630","428","558"];["r","b","g"];["L13","L15","L14"]];

% load datetimes + sats
num_events = 0;
event_fn = fullfile(direc,'event_data.txt');
fid = fopen(event_fn);
while ~feof(fid); fgetl(fid); num_events = num_events+1; end

i = 1;
times = cell(1,num_events);
sats = cell(1,num_events);
fid = fopen(event_fn);
while ~feof(fid)
    event_data = strsplit(fgetl(fid));
    times{i} = datetime([event_data{1},' ',event_data{2}]);
    sats{i} = event_data{3};
    i = i+1;
end

% load skymap
load(fullfile(direc,'dasc_data','skymap.mat'),'magnetic_footpointing')
asi_glon = magnetic_footpointing.('110km').lon;
asi_glat = magnetic_footpointing.('110km').lat;
clear('magnetic_footpointing')

for i = 1:num_events
    time = times{i};
    sat = sats{i};
    fprintf(pad('Datetime: %s - Sat: %s\n',80,'both'),time,sat)
end

fclose all;