% Description:
%   Generate a series of plots given Gemini simulation data. Plot options
%   include "closure", "conductance", "continuity", "contour", "density", "fccp",
%   "joule", "multi", or "all" for plotting every option.
%
% Example usage:
%   plot_run_2d('..\public_html\Gemini3D\maeve_4',["continuity","fccp"],300e3)
%
% Arguments:
%   direc                     gemini run directory
%   plots                     list of one or more plot name strings
%   alt_ref                   reference altitude
%   start = -1                (option) plotting start time (-1 = 0)
%   cad = -1                  (option) plotting cadence (-1 = cfg.dtout)
%   stop = -1                 (option) plotting stop time (-1 = cfg.tdur)
%   mlon_ref  = -1            (option) reference magnetic longitude (-1 = mean(MLON(:)))
%   hsv_sat  = 1e3            (option) hsv saturation magnitude
%   j_range = [-1e-6,1e-6]    (option) current standard plot range
%   n_range = [1e9,1e12]      (option) density standard plot range
%   p_range = [-3e3,3e3]      (option) potential standard plot range
%   alt_max = 400e3           (option) maximum plotting altitude
%   alt_hsv = 150e3           (option) maximum hsv plotting altitude
%   alt_cls = 120e3           (option) current closure altitude
%
% Dependencies:
%   matlab R2022a or higher
%   gemini3d
%   gemscr
%   load_conductances
%   hsv_params
%   Statistics and Machine Learning Toolbox
%   colorcet
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu

function plot_run_2d(direc,plots,alt_ref,options)
arguments
    direc (1,:) char {mustBeFolder}
    plots (1,:) string
    alt_ref (1,1) double {mustBePositive}
    options.start (1,1) double {mustBeNonempty} = -1
    options.cad (1,1) double {mustBeNonempty} = -1
    options.stop (1,1) double {mustBeNonempty} = -1
    options.mlon_ref (1,1) double {mustBeNonempty} = -1
    options.x1_range (1,2) double {mustBeNonempty} = [80e3,300e3]
    options.x3_range (1,2) double {mustBeNonempty} = [-30e3,30e3]
    options.dat (:,:) struct = struct
    %     options.j_range (1,2) double {mustBeNonempty} = [-2e-6,2e-6]
    %     options.n_range (1,2) double {mustBeNonempty} = [1e9,1e12]
    %     options.p_range (1,2) double {mustBeNonempty} = [-3e3,3e3]
end

%% assertions
plot_options = ["all","2d","conductivity","vectors"];
for plt = plots
    assert(ismember(plt,plot_options),"'"+ plt + "' is not one of the plotting options.")
end
if ismember('all',plots)
    plots = plot_options;
end
if any(direc(end)=='/\')
    direc = direc(1:end-1);
end

%% hard coded parameters
c_scl = 1e-3;       units.c = 'keV';    clm.c = 'L17';
e_scl = 1e+3;       units.e = 'mV/m';   clm.e = 'D13';
j_scl = 1e+6;       units.j = 'uA/m^2'; clm.j = 'D1A';
n_scl = 1e+0;       units.n = 'm^{-3}'; clm.n = 'L9';
s_scl = 1e+3;       units.s = 'mS/m';   clm.s = 'L18';
S_scl = 1e+0;       units.S = 'S';      clm.S = 'L18';
t_scl = 1/1.16e4;   units.t = 'eV';     clm.t = 'L3';
U_scl = 1e+3;       units.U = 'mW/m^2'; clm.U = 'L19';
v_scl = 1e-3;       units.v = 'km/s';   clm.v = 'D2';
x_scl = 1e-3;       units.x = 'km';

fts = 8*0+17; % fontsize
% ftn = 'Consolas'; % fontname (use monospaced fonts for better videos)
ftn = 'Arial';
% clb_fmt = '%+ 6.1f'; % colorbar ticklabel format
clb_fmt = '%+ 2.1f';
clb_exp = 0; % force no colorbar exponents
ctr_lc = 'k'; % contour plot linecolor
ctr_lw = 0.3; % contour plot linewidth

% j_range_hard = options.j_range;
% n_range_hard = options.n_range;
% p_range_hard = options.p_range;
% v_range_hard = [-1,1]*hsv_sat;
x1_range = options.x1_range;
x3_range = options.x3_range;
qnt = 0.99; % quantile value used to set data ranges

%% loading grid data
xg = gemini3d.read.grid(direc);
MLAT = 90-squeeze(xg.theta)*180/pi;
MLON = squeeze(xg.phi)*180/pi;
ALT = xg.alt;
x1 = xg.x1(3:end-2);
x3 = xg.x3(3:end-2);
[X1,X3] = ndgrid(x1,x3);
dx3 = xg.dx3h;
ndims = sum(xg.lx>1);
if options.mlon_ref<0
    mlon_ref = mean(MLON(:));
else
    mlon_ref = options.mlon_ref;
end

%% determining grid reference altitudes and boundaries
[~,alt_rid] = min(abs(ALT(:,1,1)-alt_ref)); % altitude reference index
[~,mlon_rid] = min(abs(MLON(1,:,1)-mlon_ref)); % magnetic longitude reference index
alt_rac = ALT(alt_rid,1,1); % altitude reference actual
mlon_rac = MLON(1,mlon_rid,1); % magnetic longitude actual

if ndims>2
    MLAT = squeeze(MLAT(:,mlon_rid,:));
end
ALT = squeeze(ALT(:,mlon_rid,:));

%% loading configuration data
cfg = gemini3d.read.config(direc);
ymd = cfg.ymd;
UTsec0 = cfg.UTsec0;
tdur = cfg.tdur;
dtout = cfg.dtout;

%% setting time boundaries
if options.start < 0
    start = 0;
else
    start = options.start;
end
if options.cad <= 0
    cad = dtout;
else
    cad = options.cad;
end
if options.stop < start
    stop = tdur;
else
    stop = options.stop;
end

%% main loop
for UTsec = UTsec0+start:cad:UTsec0+stop
    disp(pad([' UTsec = ',num2str(UTsec),' s '],80,'both','-'))

    %% loading simulation data
    time = datetime(ymd) + seconds(UTsec);
    time.Format = 'yyyyMMdd''T''HHmmss.SSS';
    if isempty(fields(options.dat))
        dat = gemini3d.read.frame(direc,'time',time);
    else
        dat = options.dat;
    end
    title_time = char(dat.time);
    filename_prefix = [char(time),'UT'];
    [~,runname] = fileparts(direc);

    %% formatting simulation data
    phi = dat.Phitop;
    [sigP,sigH,SIGP,SIGH] = tools.load_conductances(direc,time,dat,cfg,xg);

    % rescale simulation data
    j1 = dat.J1;
    j2 = dat.J2;
    j3 = dat.J3;
    ne = squeeze(dat.ne(:,mlon_rid,:));
    phi = phi(mlon_rid,:);
    SIGP = SIGP(mlon_rid,:);
    SIGH = -SIGH(mlon_rid,:); % SIGH negative?
    sigP = squeeze(sigP(:,mlon_rid,:));
    sigH = -squeeze(sigH(:,mlon_rid,:));
    Te = squeeze(dat.Te(:,mlon_rid,:));
    Ti = squeeze(dat.Ti(:,mlon_rid,:));
    v1 = dat.v1;
    v2 = dat.v2;
    v3 = dat.v3;
    if ndims>2
        j1 = squeeze(j1(:,mlon_rid,:));
        j2 = squeeze(j2(:,mlon_rid,:));
        j3 = squeeze(j3(:,mlon_rid,:));
        v1 = squeeze(v1(:,mlon_rid,:));
        v2 = squeeze(v2(:,mlon_rid,:));
        v3 = squeeze(v3(:,mlon_rid,:));
    end

    % implicit simulation data
    E3 = -gradient(phi)./dx3';

    % precipitation variables
%     [~,precip_fn,~] = fileparts(dat.filename);
%     precip = dir(fullfile(direc,'*particles',[char(precip_fn),'.*']));
%     if isempty(precip)
%         precip = dir(fullfile(direc,'*','*particles',[char(precip_fn),'.*']));
%     end
%     Q = h5read(fullfile(precip.folder,precip.name),'/Qp')/1e3; % Q defaults with mW/m^2 units
%     E0 = h5read(fullfile(precip.folder,precip.name),'/E0p'); % eV
%     Q = Q(mlon_rid,:);
%     E0 = E0(mlon_rid,:);
    %     QM = contour(squeeze(MLON(1,:,:)),squeeze(MLAT(1,:,:)),Q,1);
    %     QM = QM(:,2:end); % (mlon,mlat) contour points
    %     Q_inds = abs(QM(1,:)-mlon_rac) < 0.1*median(diff(MLON(1,:,1))); % contour line intersections
    %     Q_xlines = QM(2,Q_inds)';

    %% plotting routines
    % scale simulation data for plotting
%     E0_p = E0*c_scl; Q_p = Q*U_scl;
    E3_p = E3*e_scl;
    j1_p = j1*j_scl; j2_p = j2*j_scl; j3_p = j3*j_scl;
    ne_p = ne*n_scl;
    sigP_p = sigP*s_scl; sigH_p = sigH*s_scl;
    SIGP_p = SIGP*S_scl; SIGH_p = SIGH*S_scl;
    Te_p = Te*t_scl; Ti_p = Ti*t_scl;
    v1_p = v1*v_scl; v2_p = v2*v_scl; v3_p = v3*v_scl;

    x3_p = x3*x_scl;
    X1_p = X1*x_scl;
    X3_p = X3*x_scl;
    ALT_p = ALT*x_scl;
    alt_rac_p = round(alt_rac*x_scl);
    mlon_rac_p = round(mlon_rac);

    % set data plotting ranges
    buf = 1.05;
    j1_range_p = buf*[-1,1]*quantile(abs(j1_p(:)),qnt) + [0,1e-10];
    n_range_p = buf*[quantile(ne_p(:),1-qnt),quantile(ne_p(:),qnt)] + [0,1e-6];
    Te_range_p = buf*[quantile(Te_p(:),1-qnt),quantile(Te_p(:),qnt)] + [0,1e-6];
    Ti_range_p = buf*[quantile(Ti_p(:),1-qnt),quantile(Ti_p(:),qnt)] + [0,1e-6];
    v2_range_p = buf*[-1,1]*quantile(abs(v2_p(:)),qnt) + [0,1e-6];
    v3_range_p = buf*[-1,1]*quantile(abs(v3_p(:)),qnt) + [0,1e-6];

    %     j_range_hard_p = j_range_hard*j_scl;
    %     n_range_hard_p = n_range_hard*n_scl;
    %     p_range_hard_p = p_range_hard*p_scl;
    %     v_range_hard_p = v_range_hard*v_scl;
    x1_range_p = x1_range*x_scl;
    x3_range_p = x3_range*x_scl;

    % common plot titles and labels
    plot_title = [runname,' at ',title_time,' UT'];
    mlat_label = 'Mag. Lat.';
    north_label = ['North [',units.x,']'];
    alt_label = ['Alt. [',units.x,']'];

    reset(0)
    set(0,'defaultFigurePaperUnits','inches')
    set(0,'defaultTiledlayoutPadding','tight')
    set(0,'defaultTiledlayoutTileSpacing','tight')
    setall(0,'FontName',ftn)
    setall(0,'FontSize',fts)
    setall(0,'Multiplier',1)
    set(0,'defaultAxesFontSizeMode','manual')
    set(0,'defaultSurfaceEdgeColor','flat')
    set(0,'defaultContourLineWidth',ctr_lw)
    set(0,'defaultContourLineColor',ctr_lc)

    % all 2d plots
    if ismember('2d',plots)
        folder = '2d';
        suffix = '2d';

        figure
        set(gcf,'PaperPosition',[0,0,10,13])
        tlo = tiledlayout(5,2);
        title(tlo,[plot_title,' (mlon = ',num2str(mlon_rac_p),'°)']...
            ,'FontSize',fts...
            ,'FontName',ftn...
            ,'FontWeight'...
            ,'bold'...
            ,'Interpreter','none'...
            )

        nexttile
        pcolor(MLAT,ALT_p,log10(ne_p))
        title('Electron Density')
        %         xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.n))
        clb = colorbar;
        clb.Label.String = ['log_{10} n_e [',units.n,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(log10(n_range_p))
        ylim(x1_range_p)

        nexttile
        pcolor(MLAT,ALT_p,-j1_p)
        title('Field-Aligned Current')
        %         xlabel(mlat_label)
        %         ylabel(alt_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_{||} [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p)
        ylim(x1_range_p)

        nexttile
        pcolor(MLAT,ALT_p,j2_p)
        title('E-W Current')
        %         xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_E [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p*4)
        ylim(x1_range_p)

        nexttile
        pcolor(MLAT,ALT_p,j3_p)
        title('N-S Current')
        %         xlabel(mlat_label)
%                 ylabel(alt_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_N [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p*4)
        ylim(x1_range_p)

        nexttile
        pcolor(MLAT,ALT_p,v2_p)
        title('E-W Ion Flow')
        %         xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.v))
        clb = colorbar;
        clb.Label.String = ['v_E [',units.v,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(v2_range_p)
        ylim(x1_range_p)

        nexttile
        pcolor(MLAT,ALT_p,v3_p)
        title('N-S Ion Flow')
        %         xlabel(mlat_label)
        %         ylabel(alt_label)
        colormap(gca,colorcet(clm.v))
        clb = colorbar;
        clb.Label.String = ['v_N [',units.v,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(v3_range_p)
        ylim(x1_range_p)

        nexttile
        pcolor(MLAT,ALT_p,sigP_p)
        title('Pedersen Conductivity')
        %         xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.s))
        clb = colorbar;
        clb.Label.String = ['\sigma_P [',units.s,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        ylim(x1_range_p)

        nexttile
        pcolor(MLAT,ALT_p,sigH_p)
        title('Hall Conductivity')
        %         xlabel(mlat_label)
        %         ylabel(alt_label)
        colormap(gca,colorcet(clm.s))
        clb = colorbar;
        clb.Label.String = ['\sigma_H [',units.s,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        ylim(x1_range_p)

        nexttile
        pcolor(MLAT,ALT_p,Te_p)
        title('Electron Temperature')
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.t))
        clb = colorbar;
        clb.Label.String = ['T_e [',units.t,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(Te_range_p)
        ylim(x1_range_p)

        nexttile
        pcolor(MLAT,ALT_p,Ti_p)
        title('Ion Temperature')
        xlabel(mlat_label)
        %         ylabel(alt_label)
        colormap(gca,colorcet(clm.t))
        clb = colorbar;
        clb.Label.String = ['T_i [',units.t,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(Ti_range_p)
        ylim(x1_range_p)

        if ~exist(fullfile(direc,'plots',folder),'dir')
            mkdir(direc,fullfile('plots',folder));
        end
        filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
        disp(['Saving: ',filename])
        saveas(gcf,filename)
        close all
    end

    % conductivity plot
    if ismember('conductivity',plots)
        folder = 'conductivity';
        suffix = 'cond';

        figure
        set(gcf,'PaperPosition',[0,0,4,6])
%         title(plot_title,'FontSize',fts,'FontName',ftn,'FontWeight','bold','Interpreter','none')
%         subtitle(['alt = ',num2str(alt_rac_p),' ',units.x,', mlon = ',num2str(mlon_rac_p),'°'])

        hold on
        plot(v2_p(alt_rid,:),x3_p,'b')
        plot(-j1_p(alt_rid,:),x3_p,'r')
        plot(SIGP_p/100,x3_p,'k')
        plot(SIGH_p/100,x3_p,'--k')
        hold off
        ylabel(north_label)
        legend(...
            ['v_E [',units.v,']']...
            ,['j_{||} [',units.j,']']...
            ,['\Sigma_P [100 ',units.S,']']...
            ,['\Sigma_H [100 ',units.S,']']...
            ,'FontSize',0.5*fts,'Location','southwest','box','off')
        xlim([-1,1]*3.1)
        ylim(x3_range_p)
        grid on

        if ~exist(fullfile(direc,'plots',folder),'dir')
            mkdir(direc,fullfile('plots',folder));
        end
        filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
        disp(['Saving: ',filename])
        saveas(gcf,filename)
        close all
    end

    % vector plots
    if ismember('vectors',plots)
        folder = 'vectors';
        suffix = 'vecs';

        figure
        set(gcf,'PaperPosition',[0,0,4,6])
%         tlo = tiledlayout(1,1);
%         title(tlo,[plot_title,' (mlon = ',num2str(mlon_rac_p),'°)']...
%             ,'FontSize',fts...
%             ,'FontName',ftn...
%             ,'FontWeight'...
%             ,'bold'...
%             ,'Interpreter','none'...
%             )
        qr = 10*0+4;
        lw = 0.5;

        nexttile
        hold on
        pcolor(X3_p,X1_p,j2_p)
        quiver(...
            X3_p(1:qr/2:end,1:qr:end)...
            ,X1_p(1:qr/2:end,1:qr:end)...
            ,j3_p(1:qr/2:end,1:qr:end)...
            ,j1_p(1:qr/2:end,1:qr:end)...
            ,'-k','MarkerFaceColor','k','LineWidth',lw)
        hold off
%         title('Current Closure')
        xlabel(north_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_E [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
%         clim(j1_range_p)
        clim([-1,1]*31)
        xlim(x3_range_p)
        ylim(x1_range_p)

%         nexttile
%         hold on
%         pcolor(X3_p,X1_p,v2_p)
%         quiver(...
%             X3_p(1:qr:end,1:qr:end)...
%             ,X1_p(1:qr:end,1:qr:end)...
%             ,v3_p(1:qr:end,1:qr:end)...
%             ,v1_p(1:qr:end,1:qr:end)...
%             ,'-k','MarkerFaceColor','k','LineWidth',lw)
%         hold off
%         title('Ion Flow')
%         xlabel(north_label)
%         ylabel(alt_label)
%         colormap(gca,colorcet(clm.v))
%         clb = colorbar;
%         clb.Label.String = ['v_E [',units.v,']'];
%         clb.Ruler.TickLabelFormat = clb_fmt;
%         clb.Ruler.Exponent = clb_exp;
%         clim(v2_range_p)
%         xlim(x3_range_p)
%         ylim(x1_range_p)

        if ~exist(fullfile(direc,'plots',folder),'dir')
            mkdir(direc,fullfile('plots',folder));
        end
        filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
        disp(['Saving: ',filename])
        saveas(gcf,filename)
        close all
    end

end
    function setall(obj,property_suffix,val)
        property_suffix = char(property_suffix);
        l = length(property_suffix);
        properties = fieldnames(get(obj,'factory'));
        for n = 1:numel(properties)
            p = properties{n};
            if strcmp(p(end-l+1:end),property_suffix)
                set(obj,['default',p(8:end)],val)
            end
        end
    end
end

