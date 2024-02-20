% Description:
%   Reconstructing 2D ion flow field from distributed flow dataset
%   Uses flow data scraped from GEMINI model run, output from 'scrape.m'
%   Uses parameterized, longitudinally inclined, gaussian ridges as pseudo-basis functions
%   See reconstructor_readme.pdf for more details
%
% Arguments:
%   direc... % location of output, grid, config, and track files
%   outdir... % output files
%   orbitfn... % 'arcs_orbit_yyyymmdd.h5'
%   recon_start... % reconstruction start time, seconds from simulation start
%   recon_dur... % reconstruction duration in seconds
%   high_cad = 1            (option) cadence over time
%   droppedsats = []        (option) list of s/c to omit
%   numf = 32               (option) number of basis functions used in reconstruction
%   showplt = false         (option) show plot
%   saveplt = false         (option) save plot
%   showboundary = false    (option) show boundary plot
%   usepar = false          (option) use parallel computation
%   multirun = false        (option) when doing mutliple runs open parpool outside of function
%   maxiter = 400           (option) lsqcurvefit option value for "MaxIterations"
%   maxfuneval = 1000       (option) lsqcurvefit option value for "MaxFunctionEvaluations"
%   x2inc_lim = Inf         (option) maximum absolute potential inclination
%   scen1 = true            (option) grab arc boundary from Qp if true, else from dB data
%
% Dependencies:
%   matlab R2020b or higher
%   optimization toolbox
%   parallel computing toolbox
%   image processing toolbox
%
% Contact:
%   jules.van.irsel.GR@dartmouth.edu

function [recon_phi,recon_v2,recon_v3,P,xg] = reconstructor(direc,outdir,orbitfn,recon_start,recon_dur,options)
    arguments
        direc (1,:) char {mustBeFolder}
        outdir (1,:) char
        orbitfn (1,:) char {mustBeFile}
        recon_start (1,1) double
        recon_dur (1,1) double
        options.high_cad (1,1) int16 = 1
        options.droppedsats (1,:) int8 = []
        options.numf (1,1) int16 = 32
        options.showplt (1,1) logical = false
        options.saveplt (1,1) logical = false
        options.showboundary (1,1) logical = false
        options.usepar (1,1) logical = false
        options.multirun (1,1) logical = false
        options.maxiter (1,1) int16 = 400
        options.maxfuneval (1,1) int32 = 1000
        options.x2inc_lim (1,1) double = Inf
        options.scen1 (1,1) logical = true
    end
    global boundaryc Bmag Nm BAD %#ok<GVMIS> 
    %% basic plotting parameters
    qs = 25; % quiver plot set scale
    fz = 16; % fontsize
    ms = 1; % marker size
    lw = 1; % line width
    xl = [-800e3,800e3]; % x limits
    yl = [-300e3,300e3]; % y limits
    xt = [-800,-400,0,400,800]*1e3; % x ticks
    yt = [-300,-150,0,150,300]*1e3; % y ticks
    pc = 5; % high density plot cadence %%% CHANGE
    
    BAD = 0.0; %%% CHANGE
%     NET.addAssembly('System.Speech');
%     ssynth = System.Speech.Synthesis.SpeechSynthesizer;
%     ssynth.Volume = 100;
    
    %% grid + config metadata + loading track
    fprintf('Loading grid, configuration, and scraped data from %s...\n',direc)
    xg = gemini3d.read.grid(direc);
    cfg = gemini3d.read.config(direc);
    ymd0 = cfg.ymd;
    UTsec0 = cfg.UTsec0;
    lx2 = xg.lx(2); lx3 = xg.lx(3);
    track_fn = fullfile(direc,[orbitfn(1:end-3),'_track.mat']);
    try
        track = load(track_fn);
    catch
        error('TRACK DATA NOT FOUND IN %s',track_fn)
    end
    
    %% reconstruction parameters
    recon_stop = recon_start + recon_dur;
    full_cad = round(1/mean(gradient(track.tsat)));
    low_cad = 120;
    high_cad = options.high_cad;
    drop_sc = options.droppedsats;
    Bmag = abs(mean(xg.Bmag,'all'));
    
    %% data for comparing reconstruction
    x1alt = find(abs(xg.x1-mean(track.altsat,'all'))==min(abs(xg.x1-mean(track.altsat,'all')))); % find alt index closest to scraped alt
    dat = gemini3d.read.frame(direc, 'time', datetime([ymd0,0,0,round(UTsec0+recon_start+recon_dur/2,-1)])); % gemini output data at fraction of sim duration (pull from center of reconstruction interval)
    
    [gemini_x2,gemini_x3] = ndgrid(xg.x2(3:end-2),xg.x3(3:end-2)); % E-N distance from grid
    gemini_v2 = squeeze(dat.v2(x1alt,:,:)); % E flow at s/c alt
    gemini_v3 = squeeze(dat.v3(x1alt,:,:)); % N flow at s/c alt
    
    arcs_x2 = track.x2sat'; % E distance from s/c track
    arcs_x3 = track.x3sat'; % N distance from s/c track
    arcs_v2 = track.v2sat'; % E flow from s/c track
    arcs_v3 = track.v3sat'; % N flow from s/c track
    
    arcs_x2(:,[1 end]) = NaN; % ignore first and last datapoints...
    arcs_x3(:,[1 end]) = NaN;
    arcs_v2(:,[1 end]) = NaN;
    arcs_v3(:,[1 end]) = NaN;
    
    %% find arc boundary from Qp input
    fprintf('Fitting arc boundary...\n')
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [pth,datfn,ext] = fileparts(dat.filename);
%     dirs = dir(direc);
%     precipfn = dirs(endsWith(string({dirs.name}),'particles')).name;
%     precipfn = pth + filesep + precipfn + filesep + datfn + ext;
%     gemini_Qp = h5read(precipfn,'/Qp');
%     edges = edge(gemini_Qp,'Sobel',0.1);
%     boundary = zeros(1,lx2);
%     for ix2 = 1:lx2
%         boundary(ix2) = gemini_x3(1,find(edges(ix2,:),1,'first')); % grab southern-most detected edges
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    scen1 = options.scen1;
    if scen1 % Qp data available
        [pth,datfn,ext] = fileparts(dat.filename);
        dirs = dir(fullfile(direc,'inputs'));
        precipfn = dirs(endsWith(string({dirs.name}),'particles')).name;
        precipfn = fullfile(pth,'inputs',precipfn,datfn + ext);
        disp(precipfn)
        gemini_Qp = h5read(precipfn,'/Qp');
        edges = edge(gemini_Qp,'Sobel',0.1);    
        boundary = zeros(1,lx2);
        for ix2 = 1:lx2
            boundary(ix2) = gemini_x3(1,find(edges(ix2,:),1,'first')); % grab southern-most detected edges
        end
        boundaryf = fit(double(gemini_x2(:,1))*1e-5,double(boundary')*1e-5,'a+b*tanh(c+d*x)','Start',[0 1 0 1]); % scaling to avoid bad conditioning
    else % work with delta B data instead
%         try % TBD: streamline this. NEEDS TIME DEPENDENCE
%             mag_data = load(direc + '/mag_boundary.mat');
%         catch
%             error('MAGNETIC BOUNDARY DATA NOT FOUND')
%         end
%         p_x2 = mag_data.p_x2; % eastward distance on mag grid in meters
%         p_x3 = mag_data.p_x3; % northward distance on mag grid in meters
%         db_x2 = mag_data.db_x2; % eastward delta B 
%         db_x3 = mag_data.db_x3; % northward delta B
        
        mag_cad = 5;
        mag_t = mag_cad*round((recon_start+recon_dur/2)/mag_cad);
        osse = '2bb1'; % ASK TO CHANGE FILENAME TO OMIT THIS
        try
            p_x2 = readmatrix(direc+'\Jz_pxy_dbxy_files\Frame_t'+string(mag_t)+'px_'+osse+'.txt')*1e3; % eastward distance on mag grid in meters
            p_x3 = readmatrix(direc+'\Jz_pxy_dbxy_files\Frame_t'+string(mag_t)+'py_'+osse+'.txt')*1e3; % northward distance on mag grid in meters
            db_x2 = readmatrix(direc+'\Jz_pxy_dbxy_files\Frame_t'+string(mag_t)+'dbx_'+osse+'.txt'); % eastward delta B
            db_x3 = readmatrix(direc+'\Jz_pxy_dbxy_files\Frame_t'+string(mag_t)+'dby_'+osse+'.txt'); % northward delta B
        catch
            error('Magnetic boundary data not found for mag_t = '+string(mag_t)+'/n')
        end
            
        [~,db_max] = max(sign(db_x2).*sqrt(db_x2.^2+db_x3.^2),[],2); % location of maximum delta B magnitude per x2
        db_boundary = [p_x2(:,1)'; p_x3(1,db_max)]';
        db_boundary_x2 = db_boundary(:,1);
        db_boundary_x3 = db_boundary(:,2);
        
        edgepoints = []; % determine edge points of mag reconstruction space
        for dbi = 2:size(db_x2,1)-1
            for dbj = 2:size(db_x2,2)-1
                if ~isnan(db_x2(dbi,dbj)) % grab non-nan points only
                    if sum(isnan(db_x2(dbi,dbj+[-1 1]))) || sum(isnan(db_x2(dbi+[-1 1],dbj))) % if any adjacent point is nan, point is edge
                        edgepoints = [edgepoints; [p_x2(dbi,dbj) p_x3(dbi,dbj)]]; %#ok<AGROW> 
                    end
                end
            end
        end
        
        db_lims0 = find(diff(ismember(db_boundary,edgepoints,'rows'))==-1)+1; % find first positions where boundary points are not on edge
        db_lims1 = find(diff(ismember(db_boundary,edgepoints,'rows'))==1); % find last positions where boundary points are not yet edge
        [~,db_limi] = max(db_lims1(2:end)-db_lims0(1:end-1)); % find largest continuous interval. shift endpoints by one (always starts on edge)
        db_lim = db_lims0(db_limi)-1:db_lims1(db_limi+1)+1; % positions of largest interval without edge contact
        
        maxjump = 1e3; % smoothing parameter in meters
        db_diff = diff(db_boundary_x3);
        jumps = find(abs(db_diff) > maxjump); % positions of large jumps
        db_subtract = sum(cell2mat(arrayfun(@(x) [zeros(1,x) ones(1,size(db_x2,1)-x)]*db_diff(x),jumps,'UniformOutput',0)))'; % fill jump value to right of jump, sum for all jumps
        db_boundary_x3 = db_boundary_x3-db_subtract; % subtract the jump values from west to east
        
        db_boundary_x2 = db_boundary_x2(db_lim); % only grab boundary points within limits
        db_boundary_x3 = db_boundary_x3(db_lim);
        boundaryf = fit(db_boundary_x2*1e-5,db_boundary_x3*1e-5,'a+b*tanh(c+d*x)','Start',[0 1 0 1]); % scaling to avoid bad conditioning
    end
    
    boundaryc = coeffvalues(boundaryf).*[1e5,1e5,1,1e-5];
    arc_bound_x3 = 0*boundaryc(1)+boundaryc(2)*tanh(boundaryc(3)+boundaryc(4)*gemini_x2(:,1));

    if options.showboundary
        figure(1) % plot arc boundary
        hold on
        if scen1
            pcolor(gemini_x2/1e3,gemini_x3/1e3,gemini_Qp)
            plot(gemini_x2(:,1)/1e3,boundary/1e3,'k','LineWidth',lw)
            plot(gemini_x2(:,1)/1e3,arc_bound_x3/1e3,':r','LineWidth',lw)
        else
            pcolor(p_x2/1e3,p_x3/1e3,sign(db_x2).*sqrt(db_x2.^2+db_x3.^2))
            plot(p_x2(db_lim,1)/1e3,(db_boundary_x3-mean(db_boundary_x3))/1e3,'b','LineWidth',2*lw)
            plot(p_x2(2:end,1)/1e3,db_boundary(2:end,2)/1e3,'--k','LineWidth',2*lw)
            plot(gemini_x2(:,1)/1e3,arc_bound_x3/1e3,'--r','LineWidth',lw)
        end
        set(gcf,'units','inches','OuterPosition',[1 1 5 6])
        set(gca,'FontSize',fz)
        shading flat
%         cb0 = colorbar('southoutside');
        title('arc boundary')
        xlabel('distance east [km]')
        ylabel('distance north [km]')
%         cb0.Label.String = 'sgn(\Delta B_x)|\Delta B|';
        xlim(xl/1e3)
        ylim(yl/1e3)
        xticks(xt/1e3)
        yticks(yt/1e3)
        hold off
    end

    %% reconstruct
    istart = round(full_cad*recon_start+1);
    istop = round(full_cad*recon_stop+1);
    
    fprintf('Reconstruction setup...\n')
    bag_x2 = arcs_x2(setdiff(1 : end, drop_sc), max(istart,1) : high_cad : min(istop,end)); % 'bag of vectors'
    bag_x3 = arcs_x3(setdiff(1 : end, drop_sc), max(istart,1) : high_cad : min(istop,end));
    bag_v2 = arcs_v2(setdiff(1 : end, drop_sc), max(istart,1) : high_cad : min(istop,end));
    bag_v3 = arcs_v3(setdiff(1 : end, drop_sc), max(istart,1) : high_cad : min(istop,end));
    
    bag_x2_p = arcs_x2(setdiff(1 : 4 : 32, drop_sc), 1 : low_cad : min(istart-1,end)); % past 'bag of vectors'
    bag_x3_p = arcs_x3(setdiff(1 : 4 : 32, drop_sc), 1 : low_cad : min(istart-1,end));
    bag_v2_p = arcs_v2(setdiff(1 : 4 : 32, drop_sc), 1 : low_cad : min(istart-1,end));
    bag_v3_p = arcs_v3(setdiff(1 : 4 : 32, drop_sc), 1 : low_cad : min(istart-1,end));
    
    bag_x2_f = arcs_x2(setdiff(4 : 4 : 32, drop_sc), max(istop+1,1) : low_cad : end); % future 'bag of vectors'
    bag_x3_f = arcs_x3(setdiff(4 : 4 : 32, drop_sc), max(istop+1,1) : low_cad : end);
    bag_v2_f = arcs_v2(setdiff(4 : 4 : 32, drop_sc), max(istop+1,1) : low_cad : end);
    bag_v3_f = arcs_v3(setdiff(4 : 4 : 32, drop_sc), max(istop+1,1) : low_cad : end);
    
    xdata_d = double(reshape(cat(3,bag_x2,bag_x3),[numel(bag_x2),2])); % reshape to len(bag) x 2 array
    ydata_d = double(reshape(cat(3,bag_v2,bag_v3),[numel(bag_v2),2]));
    xdata_p = double(reshape(cat(3,bag_x2_p,bag_x3_p),[numel(bag_x2_p),2]));
    ydata_p = double(reshape(cat(3,bag_v2_p,bag_v3_p),[numel(bag_v2_p),2]));
    xdata_f = double(reshape(cat(3,bag_x2_f,bag_x3_f),[numel(bag_x2_f),2]));
    ydata_f = double(reshape(cat(3,bag_v2_f,bag_v3_f),[numel(bag_v2_f),2]));
    xdata = [xdata_p ; xdata_d ; xdata_f];
    ydata = [ydata_p ; ydata_d ; ydata_f];
    xdata(isnan(ydata(:,1)),:) = []; % remove nans from bag
    ydata(isnan(ydata(:,1)),:) = [];
%     ydata = ydata + [normrnd(0,0,size(ydata,1),1) normrnd(0,0,size(ydata,1),1)]; %%% CHANGE add measurement error
%     ydata = ydata + repmat([300 300],size(ydata,1),1); %%% CHANGE add measurement error systematic
    fprintf('Number of fitting elements: ' + string(size(ydata,1)) + '\n')

    %% optimization
    fprintf('Reconstructing flow...\n')
    if options.usepar && not(options.multirun) % start parallel pool if requested
        parpool('local');
    end
    tic
    Nm = options.numf; % number of basis functions
    P_0 = [linspace(-5,5,Nm); ones(1,Nm); zeros(1,Nm); ones(1,Nm)]; % initial parameter matrix: [x3pos (100 km), x3sig (100 km), x2inc (kV / 1e5 km), x2amp (kV)]
    lb = repmat(-[Inf Inf options.x2inc_lim Inf]',[1 Nm]); % width structure limit, TBD: somehow largest s/c sep. %%% CHANGE
    ub = repmat( [Inf Inf options.x2inc_lim Inf]',[1 Nm]);
    lsq_options = optimoptions('lsqcurvefit'...
        ,'Algorithm','levenberg-marquardt'...
        ,'UseParallel',options.usepar...
        ,'StepTolerance',1e-6...
        ,'MaxIterations',options.maxiter...
        ,'MaxFunctionEvaluations',2*numel(P_0)*options.maxfuneval...
        );
    [P,~,~,exitflag] = lsqcurvefit(@(P, xdata) F(P,xdata),P_0,xdata,ydata,lb,ub,lsq_options); % optimized parameter matrix
    recon_time = toc;
    if options.usepar && not(options.multirun)
        delete(gcp('nocreate'));
    end
    fprintf('Reconstruction time: ' + string(recon_time) + ' seconds\n')
    if exitflag <= 0
        warning('LSQCURVEFIT HAS NOT REACHED LOCAL MINIMUM')
    end

    %% creating output + plot arrays
    fprintf('Constructing output arrays...\n')
    plt_xdata = double(reshape(cat(3,gemini_x2,gemini_x3),[numel(gemini_x2),2]));
    recon_v_p = F(P,xdata_p); % past s/c locations only
    recon_v_d = F(P,xdata_d); % during s/c locations only
    recon_v_f = F(P,xdata_f); % future s/c locations only
    recon_vt = F(P,plt_xdata);
    recon_v2 = single(reshape(recon_vt(:,1),[lx2,lx3]));
    recon_v3 = single(reshape(recon_vt(:,2),[lx2,lx3]));
    recon_phi = single(phi(P,gemini_x2,gemini_x3));
    recon_phi = recon_phi - mean(recon_phi,'all');
    gemini_phi = dat.Phitop;
    gemini_phi = gemini_phi - mean(gemini_phi,'all');
    
    %% Calculating error in s/c region
    fprintf('Determining goodness of fit...\n')
    reg_bl = [arcs_x2(1 ,max(istart,1))  arcs_x3(1 ,max(istart,1))]; % 4 corners of s/c array
    reg_tl = [arcs_x2(4 ,min(istop,end)) arcs_x3(4 ,min(istop,end))];
    reg_br = [arcs_x2(29,max(istart,1))  arcs_x3(29,max(istart,1))];
    reg_tr = [arcs_x2(32,min(istop,end)) arcs_x3(32,min(istop,end))];
    regf_b = griddedInterpolant([reg_bl(1) reg_br(1)],[reg_bl(2) reg_br(2)]); % 4 boundaries of s/c array
    regf_t = griddedInterpolant([reg_tl(1) reg_tr(1)],[reg_tl(2) reg_tr(2)]);
    regf_l = griddedInterpolant([reg_bl(2) reg_tl(2)],[reg_bl(1) reg_tl(1)]);
    regf_r = griddedInterpolant([reg_br(2) reg_tr(2)],[reg_br(1) reg_tr(1)]);
    reg = false(lx2,lx3);
    for ix2 = 1:lx2
        for ix3 = 1:lx3
            x2p = gemini_x2(ix2,ix3);
            x3p = gemini_x3(ix2,ix3);
            if regf_l(x3p) < x2p && x2p < regf_r(x3p) % include coordinates left of right-most border, right of left-most, etc.
                if regf_b(x2p) < x3p && x3p < regf_t(x2p)
                    reg(ix2,ix3) = true;
                end
            end
        end
    end
    reg_d = double(reg); % for plotting purposes
    reg_d(reg_d==0) = nan;
    
    error_a_v2 = (recon_v2-gemini_v2).^2; % determine square differences in region of interest
    error_a_v3 = (recon_v3-gemini_v3).^2;
    error_p_v2 = ((recon_v2-gemini_v2)./max(abs(gemini_v2(reg)))).^2; % determine square percent error in region
    error_p_v3 = ((recon_v3-gemini_v3)./max(abs(gemini_v3(reg)))).^2;
    fprintf('Root mean square difference in v2 around s/c region is ' + string(sqrt(mean(error_a_v2(reg),'all'))) + ' m/s\n')
    fprintf('Root mean square difference in v3 around s/c region is ' + string(sqrt(mean(error_a_v3(reg),'all'))) + ' m/s\n')
    fprintf('Root mean square percent error in v2 around s/c region is ' + string(100*sqrt(mean(error_p_v2(reg),'all'))) + ' %%\n')
    fprintf('Root mean square percent error in v3 around s/c region is ' + string(100*sqrt(mean(error_p_v3(reg),'all'))) + ' %%\n')
    
    %% plotting results
    spw = 0.24;
    sph = 0.25;
    spws = (1-3*spw)/4;
    sphs = (1-2*sph)/4;
    spho = 0.07;
%     crange_phi = 1.1.*[min(gemini_phi(reg)),max(gemini_phi(reg))].*1e-3;
    crange_phi = 1.1.*[min(gemini_phi(:)),max(gemini_phi(:))].*1e-3;
%     crange_v2 = 1.1.*[min(gemini_v2(reg)),max(gemini_v2(reg))];
    crange_v2 = 1.1.*[min(gemini_v2(:)),max(gemini_v2(:))];
%     crange_v2 = 1.1*[-5000,5000]; %%% CHANGE
    crange_dv2 = [0 30];
    
    figure(2)
    set(gcf,'Visible', options.showplt)
    set(gcf,'units','inches','OuterPosition',[1 1 14 10])
    
    sp1 = subplot(2,3,1);
    hold on
    pcolor(gemini_x2/1e3, gemini_x3/1e3, gemini_phi/1e3)
    quiver(bag_x2_p/1e3,bag_x3_p/1e3,bag_v2_p*qs/1e3,bag_v3_p*qs/1e3,0,'-ob','MarkerFaceColor','b','MarkerSize',ms,'LineWidth',lw)
    quiver(bag_x2(:,1:pc:end)/1e3,bag_x3(:,1:pc:end)/1e3,bag_v2(:,1:pc:end)*qs/1e3,bag_v3(:,1:pc:end)*qs/1e3,0,'-ok','MarkerFaceColor','k','MarkerSize',ms,'LineWidth',lw)
    quiver(bag_x2_f/1e3,bag_x3_f/1e3,bag_v2_f*qs/1e3,bag_v3_f*qs/1e3,0,'-or','MarkerFaceColor','r','MarkerSize',ms,'LineWidth',lw)
    hold off
    shading flat
    cb1 = colorbar('southoutside');
    clim(crange_phi)
    set(sp1,'position',[spws 3*sphs+sph+spho spw sph])
    set(sp1,'fontsize',fz)
    title('model potential')
    xlabel('distance east [km]')
    ylabel('distance north [km]')
    cb1.Label.String = 'potential [kV]';
    xlim(xl/1e3)
    ylim(yl/1e3)
    xticks(xt/1e3)
    yticks(yt/1e3)

    sp2 = subplot(2,3,2);
    reg_x3 = regf_b(0)<gemini_x3(1,:) & gemini_x3(1,:)<regf_t(0);
    hold on
    plot(gemini_x3(lx2/2,:)/1e3,gemini_phi(lx2/2,:)/1e3,'k')
    plot(gemini_x3(lx2/2,:)/1e3,recon_phi(lx2/2,:)/1e3,'r--')
    plot(gemini_x3(lx2/2,reg_x3)/1e3,recon_phi(lx2/2,reg_x3)/1e3,'r')
    hold off
    axis([-inf inf crange_phi])
    set(sp2,'position',[2*spws+spw 3*sphs+sph+spho spw sph])
    set(sp2,'fontsize',fz)
    title('model central potential')
    xlabel('distance north [km]')
    ylabel('potential [kV]')
    xlim(yl/1e3)
    xticks(yt/1e3)

    sp3 = subplot(2,3,3);
    hold on
    pcolor(gemini_x2/1e3,gemini_x3/1e3,gemini_v2/1e3)
    if scen1; contour(gemini_x2/1e3,gemini_x3/1e3,gemini_Qp,1,'k'); end
    if ~scen1
        plot(gemini_x2(:,1)/1e3,arc_bound_x3/1e3,'k','LineWidth',lw)
        plot(gemini_x2(:,1)/1e3,boundary/1e3,'r','LineWidth',lw) %%% CHANGE
    end
    quiver(bag_x2_p/1e3,bag_x3_p/1e3,bag_v2_p*qs/1e3,bag_v3_p*qs/1e3,0,'-ob','MarkerFaceColor','b','MarkerSize',ms,'LineWidth',lw)
    quiver(bag_x2(:,1:pc:end)/1e3,bag_x3(:,1:pc:end)/1e3,bag_v2(:,1:pc:end)*qs/1e3,bag_v3(:,1:pc:end)*qs/1e3,0,'-ok','MarkerFaceColor','k','MarkerSize',ms,'LineWidth',lw)
    quiver(bag_x2_f/1e3,bag_x3_f/1e3,bag_v2_f*qs/1e3,bag_v3_f*qs/1e3,0,'-or','MarkerFaceColor','r','MarkerSize',ms,'LineWidth',lw)
    hold off
    shading flat
    cb3 = colorbar('southoutside');
    clim(crange_v2/1e3)
    set(sp3,'position',[3*spws+2*spw 3*sphs+sph+spho spw sph])
    set(sp3,'fontsize',fz)
    title('model flow east')
    xlabel('distance east [km]')
    ylabel('distance north [km]')
    cb3.Label.String = 'flow east [km/s]';
    xlim(xl/1e3)
    ylim(yl/1e3)
    xticks(xt/1e3)
    yticks(yt/1e3)

    sp4 = subplot(2,3,4);
    hold on
    pcolor(gemini_x2/1e3,gemini_x3/1e3,recon_phi/1e3)
    alpha 0.2
    pcolor(gemini_x2/1e3,gemini_x3/1e3,recon_phi.*reg_d/1e3)
    quiver(xdata_p(:,1)/1e3,xdata_p(:,2)/1e3,recon_v_p(:,1)*qs/1e3,recon_v_p(:,2)*qs/1e3,0,'-ob','MarkerFaceColor','b','MarkerSize',ms,'LineWidth',lw)
    quiver(xdata_d(1:pc:end,1)/1e3,xdata_d(1:pc:end,2)/1e3,recon_v_d(1:pc:end,1)*qs/1e3,recon_v_d(1:pc:end,2)*qs/1e3,0,'-ok','MarkerFaceColor','k','MarkerSize',ms,'LineWidth',lw)
    quiver(xdata_f(:,1)/1e3,xdata_f(:,2)/1e3,recon_v_f(:,1)*qs/1e3,recon_v_f(:,2)*qs/1e3,0,'-or','MarkerFaceColor','r','MarkerSize',ms,'LineWidth',lw)
    hold off
    shading flat
    cb4 = colorbar('southoutside');
    clim(crange_phi)
    set(sp4,'position',[spws sphs+spho spw sph])
    set(sp4,'fontsize',fz)
    title('recon. potential')
    xlabel('distance east [km]')
    ylabel('distance north [km]')
    cb4.Label.String = 'potential [kV]';
    xlim(xl/1e3)
    ylim(yl/1e3)
    xticks(xt/1e3)
    yticks(yt/1e3)

    sp5 = subplot(2,3,5);
    hold on
    pcolor(gemini_x2/1e3,gemini_x3/1e3,100.*sqrt(error_p_v2))
    alpha 0.2
    pcolor(gemini_x2/1e3,gemini_x3/1e3,100.*sqrt(error_p_v2).*reg_d)
    if scen1; contour(gemini_x2/1e3,gemini_x3/1e3,gemini_Qp,1,'k'); end
    if ~scen1
        plot(gemini_x2(:,1)/1e3,arc_bound_x3/1e3,'k','LineWidth',lw)
        plot(gemini_x2(:,1)/1e3,boundary/1e3,'r','LineWidth',lw) %%% CHANGE
    end
    hold off
    shading flat
    cb5 = colorbar('southoutside');
    clim(crange_dv2)
    set(sp5,'position',[2*spws+spw sphs+spho spw sph])
    set(sp5,'fontsize',fz)
    title('reconstruction error')
    xlabel('distance east [km]')
    ylabel('distance north [km]')
    cb5.Label.String = 'flow east error [%]';
    xlim(xl/1e3)
    ylim(yl/1e3)
    xticks(xt/1e3)
    yticks(yt/1e3)

    sp6 = subplot(2,3,6);
    hold on
    pcolor(gemini_x2/1e3,gemini_x3/1e3,recon_v2/1e3)
    alpha 0.2
    pcolor(gemini_x2/1e3,gemini_x3/1e3,recon_v2.*reg_d/1e3)
    if scen1; contour(gemini_x2/1e3,gemini_x3/1e3,gemini_Qp,1,'k'); end
    if ~scen1
        plot(gemini_x2(:,1)/1e3,arc_bound_x3/1e3,'k','LineWidth',lw)
        plot(gemini_x2(:,1)/1e3,boundary/1e3,'r','LineWidth',lw) %%% CHANGE
    end
    quiver(xdata_p(:,1)/1e3,xdata_p(:,2)/1e3,recon_v_p(:,1)*qs/1e3,recon_v_p(:,2)*qs/1e3,0,'-ob','MarkerFaceColor','b','MarkerSize',ms,'LineWidth',lw)
    quiver(xdata_d(1:pc:end,1)/1e3,xdata_d(1:pc:end,2)/1e3,recon_v_d(1:pc:end,1)*qs/1e3,recon_v_d(1:pc:end,2)*qs/1e3,0,'-ok','MarkerFaceColor','k','MarkerSize',ms,'LineWidth',lw)
    quiver(xdata_f(:,1)/1e3,xdata_f(:,2)/1e3,recon_v_f(:,1)*qs/1e3,recon_v_f(:,2)*qs/1e3,0,'-or','MarkerFaceColor','r','MarkerSize',ms,'LineWidth',lw)
    hold off
    shading flat
    cb6 = colorbar('southoutside');
    clim(crange_v2/1e3)
    set(sp6,'position',[3*spws+2*spw sphs+spho spw sph])
    set(sp6,'fontsize',fz)
    title('recon. flow east')
    xlabel('distance east [km]')
    ylabel('distance north [km]')
    cb6.Label.String = 'flow east [km/s]';
    xlim(xl/1e3)
    ylim(yl/1e3)
    xticks(xt/1e3)
    yticks(yt/1e3)
    
    if exist(fullfile(direc,outdir,'reconstructor_plots'),'dir')~=7
        mkdir(fullfile(direc,outdir,'reconstructor_plots'))
    end
    if options.saveplt
        print(fullfile(direc,outdir,'reconstructor_plots') + 'reconflow'...
            + '_' + orbitfn(12:end-3)...
            + '_' + string(full_cad)...
            + '_' + string(round(recon_start))...
            + '_' + string(round(recon_dur))...
            + '_' + string(high_cad)...
            + '_' + string(Nm)...
            + '_' + string(length(drop_sc))...
            + '_' + string(size(ydata,1))...
            + '_' + string(options.maxiter)...
            + '_' + 'P'...
            + '_' + string(round(recon_time))...
            + '_' + string(round(sqrt(mean(error_a_v2(reg),'all'))))...
            + '_' + string(round(100*sqrt(mean(error_p_v2(reg),'all'))))...
            + '_' + string(exitflag)...
            + '.png','-dpng')
    end
    
    fid = fopen(fullfile(direc,outdir,'reconstructor_report.txt'),'a');
    fprintf(fid,'%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%d\t%s\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\n',...
        datetime...
        ,direc...
        ,outdir...
        ,orbitfn...
        ,full_cad...
        ,recon_start...
        ,recon_dur...
        ,high_cad...
        ,uint16(Nm)...
        ,strjoin(['[' string(drop_sc) ']'])...
        ,uint16(size(ydata,1))...
        ,uint32(options.maxiter)...
        ,recon_time...
        ,sqrt(mean(error_a_v2(reg),'all'))...
        ,sqrt(mean(error_a_v3(reg),'all'))...
        ,100*sqrt(mean(error_p_v2(reg),'all'))...
        ,100*sqrt(mean(error_p_v3(reg),'all'))...
        ,uint8(exitflag)...
        );
    fclose(fid);

    %% functions
    % phi basis function
    function phi = phi(P,x2,x3)
        x3pos = P(1,:)*1e5;
        x3sig = P(2,:)*1e5;
        x2inc = P(3,:)*1e-4;
        x2int = P(4,:)*1e3;
        phi = 0;
        for m = 1:Nm
            b = boundaryc(1) + boundaryc(2)*tanh(boundaryc(3)+boundaryc(4)*x2(:,1)) + BAD*x2(:,1);
            phi = phi + (x2inc(m).*x2 + x2int(m)).*exp(-((x3-x3pos(m)-b)./x3sig(m)).^2);
        end
    end

    % lsqcurvefit fitting function
    function v = F(P,xdata)
        Ni = size(xdata,1); % number of vectors in bag
        x3pos = P(1,:)*1e5; % latitudinal positions [m] (P entries are near unity)
        x3sig = P(2,:)*1e5; % latitudinal widths [m]
        x2inc = P(3,:)*1e-4; % longitudenal slope of potential ridge [V/m]
        x2amp = P(4,:)*1e3; % central amplitude of potential ridge [V]
        E = zeros(Ni,2); % electric fielg
        v = zeros(Ni,2); % plasma flow field
        for i = 1:Ni % iterate through bag of vectors
            x2 = xdata(i,1); % current east position
            x3 = xdata(i,2); % current north position
            b = boundaryc(1) + boundaryc(2)*tanh(boundaryc(3)+boundaryc(4)*x2) + BAD*x2; % current boundary position
            db = boundaryc(2)*boundaryc(4)*sech(boundaryc(3)+boundaryc(4)*x2)^2 + BAD; % slope of boundary
            % caluclate the elements of -grad(phi) = -sum_m grad(phi_m) (see documentation for details)
            for m = 1:Nm % iterate through number of basis functions
                expf = exp(-((x3-x3pos(m)-b)/x3sig(m))^2);
                E(i,1) = E(i,1) + (-x2inc(m)-(2/x3sig(m)^2)*(x2inc(m)*x2+x2amp(m))*(x3-x3pos(m)-b)*db)*expf;
                E(i,2) = E(i,2) + (2/x3sig(m)^2)*(x2inc(m)*x2+x2amp(m))*(x3-x3pos(m)-b)*expf;
            end
            v(:,1) = -E(:,2)./Bmag; % v = ExB/B^2
            v(:,2) =  E(:,1)./Bmag;
        end
    end
end
