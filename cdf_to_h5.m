direc = fullfile('data','swop','swarm_data');
cdfs = {dir(fullfile(direc,'*.cdf')).name};

for cdf = cdfs
    path = fullfile(direc,cdf{1});
    path_h5 = [path(1:end-3),'h5'];
    [all_data, info] = cdfread(path,'DatetimeType','datetime'); 
    vars = info.Variables(:,1);
    units = info.VariableAttributes.UNITS(:,2);
    types = info.VariableAttributes.Type(:,2);
    descs = info.VariableAttributes.CATDESC(:,2);

    for i = 1:length(vars)
        data = [all_data{:,i}];
        var = ['/',vars{i}];
        unit = units{i};
        desc = descs{i};
        type = types{i};
        type = lower(type(5:end));
        if strcmp(type,'float')
            type = 'double';
        elseif strcmp(type(2:end-1),"int")
            type = [type(1:end-1),'32'];
        elseif strcmp(type,'epoch')
            type = 'double';
            data = posixtime(data);
            unit = 's';
            desc = 'Unix time.';
        end
        try
            h5create(path_h5,var,size(data),'Datatype',type)
            h5write(path_h5,var,data)
            h5writeatt(path_h5,var,'Units',unit)
            h5writeatt(path_h5,var,'Description',desc)
        catch
            warning('%s in %s already exists',var(2:end),path_h5)
        end

        if strcmp(var,'/Latitude'); lat_gc = data; end
        if strcmp(var,'/Radius'); rad = data; end
    end
    
    global_atts = fieldnames(info.GlobalAttributes)';
    for att = global_atts
        h5writeatt(path_h5,'/',att{1},info.GlobalAttributes.(att{1}))
    end

    [lat_gd,alt_gd] = geoc2geod(lat_gc,rad);
    h5create(path_h5,'/GeodeticLatitude',size(lat_gd),'Datatype','double')
    h5write(path_h5,'/GeodeticLatitude',lat_gd)
    h5writeatt(path_h5,'/GeodeticLatitude','Units','degrees')
    h5writeatt(path_h5,'/GeodeticLatitude','Description','Geodetic Latitude.')
    h5create(path_h5,'/GeodeticAltitude',size(alt_gd),'Datatype','double')
    h5write(path_h5,'/GeodeticAltitude',alt_gd)
    h5writeatt(path_h5,'/GeodeticAltitude','Units','m')
    h5writeatt(path_h5,'/GeodeticAltitude','Description','Geodetic Altitude.')
end