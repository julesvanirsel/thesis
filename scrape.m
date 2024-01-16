function track = scrape(direc)
    arguments
        direc (1,:) char {mustBeFolder}
    end
    
    %% grid + config metadata
    fprintf('Loading grid and configuration data from %s\n',direc)
    cfg = gemini3d.read.config(direc);
    UTsec0 = cfg.UTsec0;
    tdur = cfg.tdur;

    %% satilite orbit data
    orbitfn = 'arcs_orbit_20210211.h5';
    fprintf('Loading orbit data from %s\n',orbitfn)
    dates = h5read(orbitfn,'/ModJDate');
    glatsat = h5read(orbitfn,'/Lat');
    glonsat = h5read(orbitfn,'/Lon');
    altsat = h5read(orbitfn,'/Alt');
    tsat = (dates-dates(1)).*(24*60*60); % seconds from start of day

    %% limit orbit times to simulation times
    inds = round(tsat) >= UTsec0 & round(tsat) <= UTsec0+tdur;
    glatsat = glatsat(inds,:);
    glonsat = glonsat(inds,:);
    altsat = altsat(inds,:);
    tsat = tsat(inds)-UTsec0; % seconds from start of sim

    %% interpolate orbit data to new cadence
    fprintf('Interpolating orbit data to new cadence\n')
    cad = 30; % approximate cadence in Hz
    tsatnew = linspace(min(tsat),max(tsat),round(range(tsat)*cad));
    lsat = length(glatsat);
    glatsatnew = zeros(length(tsatnew),lsat);
    glonsatnew = zeros(length(tsatnew),lsat);
    altsatnew = zeros(length(tsatnew),lsat);
    for isat = 1:lsat
        glatf = griddedInterpolant(tsat,glatsat(:,isat));
        glonf = griddedInterpolant(tsat,wrapTo360(glonsat(:,isat))); % 0-360 deg as quick fix for making interpolation continuous 
        altf = griddedInterpolant(tsat,altsat(:,isat));
        glatsatnew(:,isat) = glatf(tsatnew);
        glonsatnew(:,isat) = wrapTo180(glonf(tsatnew)); % return to -180 - 180 deg
        altsatnew(:,isat) = altf(tsatnew);
    end

    %% scrape data
    tic
    fprintf('Scraping orbits through model space\n')
    track = gemscr.postprocess.virtual_spacecraft_cartesian(direc,glonsatnew,glatsatnew,altsatnew,tsatnew);
    scrape_toc = toc;
    fprintf('Scraping time: ' + string(scrape_toc) + ' seconds.\n')
    save(fullfile(direc,[orbitfn(1:end-3),'_track.mat']),'-struct','track')
end
