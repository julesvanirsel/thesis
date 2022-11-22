% Description:
%   Generate a series of plots given Gemini simulation data. Plot options
%   include "closure", "conductance", "continuity", "contour", "density", "fccp",
%   "joule", "multi", or "all" for plotting every option.
%
% Example usage:
%   plot_run('..\public_html\Gemini3D\maeve_4',["continuity","fccp"],300e3)
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

function plot_run(direc,plots,alt_ref,options)
arguments
    direc (1,:) char {mustBeFolder}
    plots (1,:) string
    alt_ref (1,1) double {mustBePositive}
    options.start (1,1) double {mustBeNonempty} = -1
    options.cad (1,1) double {mustBeNonempty} = -1
    options.stop (1,1) double {mustBeNonempty} = -1
    options.mlon_ref (1,1) double {mustBeNonempty} = -1
    options.hsv_sat (1,1) double {mustBeNonempty} = 1e3
    options.j_range (1,2) double {mustBeNonempty} = [-2e-6,2e-6]
    options.n_range (1,2) double {mustBeNonempty} = [1e9,1e12]
    options.p_range (1,2) double {mustBeNonempty} = [-3e3,3e3]
    options.alt_max (1,1) double {mustBeNonempty} = 400e3
    options.alt_hsv (1,1) double {mustBeNonempty} = 150e3
    options.alt_cls (1,1) double {mustBeNonempty} = 120e3
end

%% assertions
plot_options = ["all","closure","conductance","continuity","contour","density","fccp","joule","multi","temp"];
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
c_scl = 1e-3; units.c = 'keV';    clm.c = 'L17';
e_scl = 1e+3; units.e = 'mV/m';   clm.e = 'D13';
j_scl = 1e+6; units.j = 'uA/m^2'; clm.j = 'D1A';
n_scl = 1e+0; units.n = 'm^{-3}'; clm.n = 'L9';
p_scl = 1e-3; units.p = 'kV';     clm.p = 'D10';
s_scl = 1e+0; units.s = 'S';      clm.s = 'L18';
t_scl = 1e-3; units.t = 'kK';     clm.t = 'L3';
u_scl = 1e+6; units.u = 'uW/m^3'; clm.u = 'L19';
U_scl = 1e+3; units.U = 'mW/m^2'; clm.U = 'L19';
v_scl = 1e-3; units.v = 'km/s';   clm.v = 'D2';
x_scl = 1e-3; units.x = 'km';

fts = 8; % fontsize
ftn = 'Consolas'; % fontname (use monospaced fonts for better videos)
clb_fmt = '%+ 6.1f'; % colorbar ticklabel format
clb_exp = 0; % force no colorbar exponents
ctr_lc = 'k'; % contour plot linecolor
ctr_lw = 0.3; % contour plot linewidth

hsv_sat = options.hsv_sat;
j_range_hard = options.j_range;
n_range_hard = options.n_range;
p_range_hard = options.p_range;
v_range_hard = [-1,1]*hsv_sat;
alt_max = options.alt_max;
alt_hsv = options.alt_hsv;
alt_cls = options.alt_cls;
qnt = 0.99; % quantile value used to set data ranges

%% loading grid data
xg = gemini3d.read.grid(direc);
MLAT = 90-squeeze(xg.theta)*180/pi;
MLON = squeeze(xg.phi)*180/pi;
ALT = xg.alt;
x2 = xg.x2(3:end-2);
x3 = xg.x3(3:end-2);
[X2,X3] = ndgrid(x2,x3);
lx2 = xg.lx(2); lx3 = xg.lx(3);
dx1 = xg.dx1h;
if options.mlon_ref<0
    mlon_ref = mean(MLON(:));
else
    mlon_ref = options.mlon_ref;
end

%% determining grid reference altitudes and boundaries
lb1 = 1;
[~,ub1] = min(abs(ALT(:,1,1)-alt_max));
lb2 = 3;
ub2 = lx2+1-lb2;
lb3 = 3;
ub3 = lx3+1-lb3;
MLAT = MLAT(lb1:ub1,lb2:ub2,lb3:ub3);
MLON = MLON(lb1:ub1,lb2:ub2,lb3:ub3);
ALT = ALT(lb1:ub1,lb2:ub2,lb3:ub3);
X2 = X2(lb2:ub2,lb3:ub3);
X3 = X3(lb2:ub2,lb3:ub3);
dx1 = dx1(lb1:ub1);
[~,alt_rid] = min(abs(ALT(:,1,1)-alt_ref)); % altitude reference index
[~,alt_cid] = min(abs(ALT(:,1,1)-alt_cls)); % altitude closure index
[~,mlon_rid] = min(abs(MLON(1,:,1)-mlon_ref)); % magnetic longitude reference index
alt_rac = ALT(alt_rid,1,1); % altitude reference actual
alt_cac = ALT(alt_cid,1,1); % altitude closure actual
mlon_rac = MLON(1,mlon_rid,1); % magnetic longitude actual

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
    dat = gemini3d.read.frame(direc,'time',time);
    title_time = char(dat.time);
    filename_prefix = [char(time),'UT'];
    [~,runname] = fileparts(direc);

    %% formatting simulation data
    phi = dat.Phitop;
    [E1,E2,E3] = gemscr.postprocess.pot2field(xg,phi);
    [jP_3,jH_3,~] = gemscr.postprocess.current_decompose(xg,dat);
    [sigP,sigH,SIGP,SIGH] = load_conductances(direc,time,dat,cfg,xg);

    % rescale simulation data
    E1 = E1(lb1:ub1,lb2:ub2,lb3:ub3);
    E2 = E2(lb1:ub1,lb2:ub2,lb3:ub3);
    E3 = E3(lb1:ub1,lb2:ub2,lb3:ub3);
    j1 = dat.J1(lb1:ub1,lb2:ub2,lb3:ub3);
    j2 = dat.J2(lb1:ub1,lb2:ub2,lb3:ub3);
    j3 = dat.J3(lb1:ub1,lb2:ub2,lb3:ub3);
    jP_3 = jP_3(lb1:ub1,lb2:ub2,lb3:ub3,:);
    jH_3 = jH_3(lb1:ub1,lb2:ub2,lb3:ub3,:);
    jP_3 = sqrt(sum(jP_3.^2,4)); % get magnitudes
    jH_3 = sqrt(sum(jH_3.^2,4));
    ne = dat.ne(lb1:ub1,lb2:ub2,lb3:ub3);
    phi = phi(lb2:ub2,lb3:ub3);
    SIGP = SIGP(lb2:ub2,lb3:ub3);
    SIGH = -SIGH(lb2:ub2,lb3:ub3); % SIGH negative?
    sigP = sigP(lb1:ub1,lb2:ub2,lb3:ub3);
    sigH = -sigH(lb1:ub1,lb2:ub2,lb3:ub3);
    Te = dat.Te(lb1:ub1,lb2:ub2,lb3:ub3);
    Ti = dat.Ti(lb1:ub1,lb2:ub2,lb3:ub3);
    v2 = dat.v2(lb1:ub1,lb2:ub2,lb3:ub3);
    v3 = dat.v3(lb1:ub1,lb2:ub2,lb3:ub3);

    % implicit simulation data
    jP = (j2.*E2+j3.*E3)./sqrt(E2.^2+E3.^2); % E_perp direction
    jH = (j2.*E3-j3.*E2)./sqrt(E2.^2+E3.^2); % b x E_perp direction
    jP_2 = sigP.*sqrt(E2.^2+E3.^2);
    jH_2 = sigH.*sqrt(E2.^2+E3.^2);
    [~,dX2] = gradient(X2);
    [dX3,~] = gradient(X3);
    [~,dE22] = gradient(squeeze(E2(1,:,:)));
    [dE33,~] = gradient(squeeze(E3(1,:,:)));
    [dSIGP3,dSIGP2] = gradient(SIGP);
    [dSIGH3,dSIGH2] = gradient(SIGH);
    jA = SIGP.*(dE22./dX2+dE33./dX3); % SIGP*div_perp(E)
    jB = (dSIGP2.*squeeze(E2(1,:,:))./dX2 + dSIGP3.*squeeze(E3(1,:,:))./dX3); % grad(SIGP).E_perp
    jC = (dSIGH2.*squeeze(E3(1,:,:))./dX2 - dSIGH3.*squeeze(E2(1,:,:))./dX3); % grad(SIGH).bxE_perp
    joule = (j1.*E1 + j2.*E2 + j3.*E3);
    joule_int = squeeze(sum(joule.*dx1,1));
    joule_int_2 = (E1.*E1 + E2.*E2 + E3.*E3);
    joule_int_2 = SIGP.*squeeze(joule_int_2(1,:,:));

    % hsv plot variables
    [hsv_map_clb,hsv_mlon,hsv_mlon_map,hsv_alt,hsv_alt_map] =...
        hsv_params(v2,v3,MLAT,MLON,ALT,alt_ref,mlon_ref,hsv_sat);

    % precipitation variables
    [~,precip_fn,~] = fileparts(dat.filename);
    precip = dir(fullfile(direc,'*particles',[char(precip_fn),'.*']));
    if isempty(precip)
        precip = dir(fullfile(direc,'*','*particles',[char(precip_fn),'.*']));
    end
    Q = h5read(fullfile(precip.folder,precip.name),'/Qp')/1e3; % Q defaults with mW/m^2 units
    E0 = h5read(fullfile(precip.folder,precip.name),'/E0p'); % eV
    Q = Q(lb2:ub2,lb3:ub3);
    E0 = E0(lb2:ub2,lb3:ub3);
    QM = contour(squeeze(MLON(1,:,:)),squeeze(MLAT(1,:,:)),Q,1);
    QM = QM(:,2:end); % (mlon,mlat) contour points
    Q_inds = abs(QM(1,:)-mlon_rac) < 0.1*median(diff(MLON(1,:,1))); % contour line intersections
    Q_xlines = QM(2,Q_inds)';

    %% plotting routines
    % scale simulation data for plotting
    E0_p = E0*c_scl; Q_p = Q*U_scl;
    E2_p = E2*e_scl; E3_p = E3*e_scl;
    j1_p = j1*j_scl; jP_p = jP*j_scl; jH_p = jH*j_scl;
    jP_2_p = jP_2*j_scl; jH_2_p = jH_2*j_scl; jP_3_p = jP_3*j_scl; jH_3_p = jH_3*j_scl;
    jA_p = jA*j_scl; jB_p = jB*j_scl; jC_p = jC*j_scl;
    ne_p = ne*n_scl;
    phi_p = phi*p_scl;
    SIGP_p = SIGP*s_scl; SIGH_p = SIGH*s_scl;
    Te_p = Te*t_scl; Ti_p = Ti*t_scl;
    joule_p = joule*u_scl;
    joule_int_p = joule_int*U_scl; joule_int_2_p = joule_int_2*U_scl;
    v2_p = v2*v_scl; v3_p = v3*v_scl;
    hsv_sat_p = hsv_sat*v_scl;

    ALT_p = ALT*x_scl;
    alt_hsv_p = alt_hsv*x_scl;
    alt_rac_p = round(alt_rac*x_scl);
    alt_cac_p = round(alt_cac*x_scl);
    mlon_rac_p = round(mlon_rac);

    % set data plotting ranges
    buf = 1.05;
    e_range_p = buf*[-1,1]*quantile(abs([E2_p(:);E3_p(:)]),qnt);
    j1_range_p = buf*[-1,1]*quantile(abs(j1_p(:)),qnt);
    jP_range_p = buf*[-1,1]*quantile(abs(jP_p(:)),qnt);
    jH_range_p = buf*[-1,1]*quantile(abs(jH_p(:)),qnt);
    n_range_p = buf*[quantile(ne_p(:),1-qnt),quantile(ne_p(:),qnt)];
    Te_range_p = buf*[quantile(Te_p(:),1-qnt),quantile(Te_p(:),qnt)];
    Ti_range_p = buf*[quantile(Ti_p(:),1-qnt),quantile(Ti_p(:),qnt)];
    v2_range_p = buf*[-1,1]*quantile(abs(v2_p(:)),qnt);
    v3_range_p = buf*[-1,1]*quantile(abs(v3_p(:)),qnt);
    u_range_p = buf*[quantile(joule_p(:),1-qnt),quantile(joule_p(:),qnt)];
    u_range_int_p = buf*[quantile(joule_int_p(:),1-qnt),quantile(joule_int_p(:),qnt)];
    p_range_p = buf*[-1,1]*quantile(abs(phi_p(:)),qnt);

    j_range_hard_p = j_range_hard*j_scl;
    n_range_hard_p = n_range_hard*n_scl;
    p_range_hard_p = p_range_hard*p_scl;
    v_range_hard_p = v_range_hard*v_scl;

    % common plot titles and labels
    plot_title = [runname,' at ',title_time,' UT'];
    mlon_label = 'Mag. Lon.';
    mlat_label = 'Mag. Lat.';
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

    % closure current plot
    if ismember('closure',plots)
        folder = 'closure';
        suffix = 'clos';

        figure
        set(gcf,'PaperPosition',[0,0,6.5,4.5])
        tlo = tiledlayout(3,2);
        title(tlo,plot_title,'FontSize',fts,'FontName',ftn,'FontWeight','bold','Interpreter','none')

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(jP_p(:,mlon_rid,:)))
        title(['Pedersen Current 1 (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_P [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(jP_range_p)
        ylim([min(ALT_p(:)),alt_hsv_p])

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(jH_p(:,mlon_rid,:)))
        title(['Hall Current 1 (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_P [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(jH_range_p)
        ylim([min(ALT_p(:)),alt_hsv_p])

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(jP_2_p(:,mlon_rid,:)))
        title(['Pedersen Current 2 (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_P [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(jP_range_p)
        ylim([min(ALT_p(:)),alt_hsv_p])

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(jH_2_p(:,mlon_rid,:)))
        title(['Hall Current 2 (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_P [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(jH_range_p)
        ylim([min(ALT_p(:)),alt_hsv_p])

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(jP_3_p(:,mlon_rid,:)))
        title(['Pedersen Current 3 (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_P [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(jP_range_p)
        ylim([min(ALT_p(:)),alt_hsv_p])

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(jH_3_p(:,mlon_rid,:)))
        title(['Hall Current 3 (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_P [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(jH_range_p)
        ylim([min(ALT_p(:)),alt_hsv_p])

        if ~exist(fullfile(direc,'plots',folder),'dir')
            mkdir(direc,fullfile('plots',folder));
        end
        filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
        disp(['Saving: ',filename])
        saveas(gcf,filename)
        close all
    end

    % conductance plot
    if ismember('conductance',plots)
        folder = 'conductance';
        suffix = 'cond';

        figure
        set(gcf,'PaperPosition',[0,0,6.5,1.5])
        tlo = tiledlayout(1,2);
        title(tlo,plot_title,'FontSize',fts,'FontName',ftn,'FontWeight','bold','Interpreter','none')

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),abs(SIGP_p))
        title('Pedersen Conductance')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.s))
        clb = colorbar;
        clb.Label.String = ['\int_{||} |\sigma_P| [',units.s,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),abs(SIGH_p))
        title('Hall Conductance')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.s))
        clb = colorbar;
        clb.Label.String = ['\int_{||} |\sigma_H| [',units.s,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;

        if ~exist(fullfile(direc,'plots',folder),'dir')
            mkdir(direc,fullfile('plots',folder));
        end
        filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
        disp(['Saving: ',filename])
        saveas(gcf,filename)
        close all
    end

    % continuity plot
    if ismember('continuity',plots)
        folder = 'continuity';
        suffix = 'cont';

        figure
        set(gcf,'PaperPosition',[0,0,6.5,4.5])
        tlo = tiledlayout(3,2);
        title(tlo,plot_title,'FontSize',fts,'FontName',ftn,'FontWeight','bold','Interpreter','none')

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),jA_p)
        title('div(E_\perp) FAC Term')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['\Sigma_P \nabla\cdot E_\perp [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p)

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),jA_p+jB_p+jC_p)
        title('Sum of FAC Terms')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_{||} [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p)

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),jB_p)
        title('grad(\Sigma_P) FAC Term')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['\nabla\Sigma_P\cdot E_\perp [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p)

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),squeeze(-j1_p(end,:,:)))
        title('Model FAC')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_{||} [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p)

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),jC_p)
        title('grad(\Sigma_H) FAC Term')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['\nabla\Sigma_H\cdot b\times E_\perp [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p)

        nexttile
        va = squeeze(-j1_p(end,:,:));
        ve = jA_p+jB_p+jC_p;
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),va-ve)
        title('FAC Difference')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['\Delta j_{||} [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p)

        if ~exist(fullfile(direc,'plots',folder),'dir')
            mkdir(direc,fullfile('plots',folder));
        end
        filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
        disp(['Saving: ',filename])
        saveas(gcf,filename)
        close all
    end

    % contour plot
    if ismember('contour',plots)
        for isauto = [0,1]
            if isauto
                folder = 'contour-auto';
                title_suffix = ' (autocolored)';
            else
                folder = 'contour-standard';
                title_suffix = '';
            end
            suffix = 'cntr';

            figure
            set(gcf,'PaperPosition',[0,0,9,3.5])
            tlo = tiledlayout(3,3);
            title(tlo,[plot_title,title_suffix],'FontSize',fts,'FontName',ftn,'FontWeight','bold','Interpreter','none')

            nexttile(1,[2,1])
            hold on
            pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),hsv_alt);
            contour(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),Q_p,1)
            title(['Ion Flow (',num2str(alt_rac_p),' ',units.x,')'])
            xlabel(mlon_label)
            ylabel(mlat_label)
            colormap(gca,hsv_alt_map)
            clb = colorbar;
            colormap(clb,hsv_map_clb)
            clb.Limits = [0,1];
            clb.Ticks = [0,1/4,1/2,3/4,1];
            clb.TickLabels = {'W','S','E','N','W'};
            clb.Label.String = ['Sat. at ',num2str(hsv_sat_p),' ',units.v];
            clim([0,1])
            xline(mlon_rac_p,'r--')

            nexttile(2,[2,1])
            hold on
            pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),phi_p)
            contour(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),Q_p,1)
            title('Electric Potential')
            xlabel(mlon_label)
            ylabel(mlat_label)
            colormap(gca,colorcet(clm.p))
            clb = colorbar;
            clb.Label.String = ['\phi_{top} [',units.p,']'];
            clb.Ruler.TickLabelFormat = clb_fmt;
            clb.Ruler.Exponent = clb_exp;
            if isauto
                clim(p_range_p)
            else
                clim(p_range_hard_p)
            end

            nexttile(3,[2,1])
            hold on
            pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(-j1_p(alt_rid,:,:)))
            contour(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),Q_p,1)
            title(['FAC (',num2str(alt_rac_p),' ',units.x,')'])
            xlabel(mlon_label)
            ylabel(mlat_label)
            colormap(gca,colorcet(clm.j))
            clb = colorbar;
            clb.Label.String = ['J_{||} [',units.j,']'];
            clb.Ruler.TickLabelFormat = clb_fmt;
            clb.Ruler.Exponent = clb_exp;
            if isauto
                clim(j1_range_p)
            else
                clim(j_range_hard_p)
            end
            xline(mlon_rac_p,'r--')

            nexttile
            pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),hsv_mlon)
            title(['Ion Flow (',num2str(mlon_rac_p),'°)'])
            xlabel(mlat_label)
            ylabel(alt_label)
            colormap(gca,hsv_mlon_map)
            clb = colorbar;
            colormap(clb,hsv_map_clb)
            clb.Limits = [0,1];
            clb.Ticks = [0,1/4,1/2,3/4,1];
            clb.TickLabels = {'W','S','E','N','W'};
            clb.Label.String = ['Sat. at ',num2str(hsv_sat_p),' ',units.v];
            clim([0,1])
            ylim([min(ALT_p(:)),alt_hsv_p])
            if ~isempty(Q_xlines)
                xline(Q_xlines,'LineWidth',ctr_lw)
            end
            yline(alt_rac_p,'r--')

            nexttile
            pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(log10(ne_p(:,mlon_rid,:))))
            title(['Density (',num2str(mlon_rac_p),'°)'])
            xlabel(mlat_label)
            ylabel(alt_label)
            colormap(gca,colorcet(clm.n))
            clb = colorbar;
            clb.Label.String = ['log_{10} n_e [',units.n,']'];
            clb.Ruler.TickLabelFormat = clb_fmt;
            clb.Ruler.Exponent = clb_exp;
            if isauto
                clim(log10(n_range_p))
            else
                clim(log10(n_range_hard_p))
            end
            if ~isempty(Q_xlines)
                xline(Q_xlines,'LineWidth',ctr_lw)
            end

            nexttile
            pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(-j1_p(:,mlon_rid,:)))
            title(['FAC (',num2str(mlon_rac_p),'°)'])
            xlabel(mlat_label)
            ylabel(alt_label)
            colormap(gca,colorcet(clm.j))
            clb = colorbar;
            clb.Label.String = ['j_{||} [',units.j,']'];
            clb.Ruler.TickLabelFormat = clb_fmt;
            clb.Ruler.Exponent = clb_exp;
            if isauto
                clim(j1_range_p)
            else
                clim(j_range_hard_p)
            end
            if ~isempty(Q_xlines)
                xline(Q_xlines,'LineWidth',ctr_lw)
            end
            yline(alt_rac_p,'r--')

            if ~exist(fullfile(direc,'plots',folder),'dir')
                mkdir(direc,fullfile('plots',folder));
            end
            filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
            disp(['Saving: ',filename])
            saveas(gcf,filename)
            close all
        end
    end

    % density plot
    if ismember('density',plots)
        folder = 'density';
        suffix = 'dens';

        figure
        set(gcf,'PaperPosition',[0,0,6.5,3])
        tlo = tiledlayout(2,2);
        title(tlo,plot_title,'FontSize',fts,'FontName',ftn,'FontWeight','bold','Interpreter','none')

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(log10(ne_p(alt_rid,:,:))))
        title(['Density (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.n))
        clb = colorbar;
        clb.Label.String = ['log_{10} n_e [',units.n,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(log10(n_range_p))
        xline(mlon_rac_p,'r--')

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(log10(ne_p(:,mlon_rid,:))))
        title(['Density (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.n))
        clb = colorbar;
        clb.Label.String = ['log_{10} n_e [',units.n,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(log10(n_range_p))
        yline(alt_rac_p,'r--')

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(-j1_p(alt_rid,:,:)))
        title(['FAC (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['j_{||} [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p)

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(v2_p(alt_rid,:,:)))
        title(['Ion Flow (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.v))
        clb = colorbar;
        clb.Label.String = ['v_E [',units.v,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(v2_range_p)

        if ~exist(fullfile(direc,'plots',folder),'dir')
            mkdir(direc,fullfile('plots',folder));
        end
        filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
        disp(['Saving: ',filename])
        saveas(gcf,filename)
        close all
    end

    % flow-current-conductance-precipitation plot
    if ismember('fccp',plots)
        folder = 'fccp';
        suffix = 'fccp';

        figure
        set(gcf,'PaperPosition',[0,0,6.5,6])
        tlo = tiledlayout(4,2);
        title(tlo,plot_title,'FontSize',fts,'FontName',ftn,'FontWeight','bold','Interpreter','none')

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(v2_p(alt_rid,:,:)))
        title(['E-W Ion Flow (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.v))
        clb = colorbar;
        clb.Label.String = ['v_E [',units.v,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(v2_range_p)

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(v3_p(alt_rid,:,:)))
        title(['N-S Ion Flow (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.v))
        clb = colorbar;
        clb.Label.String = ['v_N [',units.v,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(v3_range_p)

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),hsv_alt);
        title(['Ion Flow (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,hsv_alt_map)
        clb = colorbar;
        colormap(clb,hsv_map_clb)
        clb.Limits = [0,1];
        clb.Ticks = [0,1/4,1/2,3/4,1];
        clb.TickLabels = {'W','S','E','N','W'};
        clb.Label.String = ['Sat. at ',num2str(hsv_sat_p),' ',units.v];
        clim([0,1])

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(-j1_p(alt_rid,:,:)))
        title(['FAC (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['J_{||} [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p)

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),abs(SIGP_p))
        title('Pedersen Conductance')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.s))
        clb = colorbar;
        clb.Label.String = ['\int_{||} |\sigma_P| [',units.s,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),abs(SIGH_p))
        title('Hall Conductance')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.s))
        clb = colorbar;
        clb.Label.String = ['\int_{||} |\sigma_H| [',units.s,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),Q_p)
        title('Precip. Energy Flux')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.U))
        clb = colorbar;
        clb.Label.String = ['Q [',units.U,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),E0_p)
        title('Charactaristic Energy')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.c))
        clb = colorbar;
        clb.Label.String = ['E_0 [',units.c,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;

        if ~exist(fullfile(direc,'plots',folder),'dir')
            mkdir(direc,fullfile('plots',folder));
        end
        filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
        disp(['Saving: ',filename])
        saveas(gcf,filename)
        close all
    end

    % joule heating plot
    if ismember('joule',plots)
        folder = 'jouleheating';
        suffix = 'joul';

        figure
        set(gcf,'PaperPosition',[0,0,6.5,6])
        tlo = tiledlayout(4,2);
        title(tlo,plot_title,'FontSize',fts,'FontName',ftn,'FontWeight','bold','Interpreter','none')

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(E2_p(alt_rid,:,:)))
        title(['E-W Elec. Field (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.e))
        clb = colorbar;
        clb.Label.String = ['E_E [',units.e,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(e_range_p)

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(E3_p(alt_rid,:,:)))
        title(['N-S Elec. Field (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.e))
        clb = colorbar;
        clb.Label.String = ['E_N [',units.e,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(e_range_p)

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(-j1_p(alt_rid,:,:)))
        title(['FAC (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['J_{||} [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p)
        xline(mlon_rac_p,'r--')

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(-j1_p(:,mlon_rid,:)))
        title(['FAC (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.j))
        clb = colorbar;
        clb.Label.String = ['J_{||} [',units.j,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(j1_range_p)
        yline(alt_rac_p,'r--')

        nexttile
        pcolor(squeeze(MLON(alt_cid,:,:)),squeeze(MLAT(alt_cid,:,:)),squeeze(joule_p(alt_cid,:,:)))
        title(['Joule Heating (',num2str(alt_cac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.u))
        clb = colorbar;
        clb.Label.String = ['J\cdot E [',units.u,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(u_range_p)
        xline(mlon_rac_p,'r--')

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(joule_p(:,mlon_rid,:)))
        title(['Joule Heating (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.u))
        clb = colorbar;
        clb.Label.String = ['J\cdot E [',units.u,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(u_range_p)
        yline(alt_cac_p,'r--')

        nexttile
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),joule_int_p)
        title('Integrated Joule Heating')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.U))
        clb = colorbar;
        clb.Label.String = ['\int_{||} J\cdot E [',units.U,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(u_range_int_p)
        xline(mlon_rac_p,'r--')

        nexttile
        va = joule_int_p;
        ve = joule_int_2_p;
        pcolor(squeeze(MLON(end,:,:)),squeeze(MLAT(end,:,:)),100*(va./ve-1))
        title('Percent Error')
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet('D1'))
        clb = colorbar;
        clb.Label.String = '\int_{||} J\cdot E/(\Sigma_P E^2)-1 [%]';
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim([-1,1]*5)
        xline(mlon_rac_p,'r--')

        if ~exist(fullfile(direc,'plots',folder),'dir')
            mkdir(direc,fullfile('plots',folder));
        end
        filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
        disp(['Saving: ',filename])
        saveas(gcf,filename)
        close all
    end

    % multipanel plot
    if ismember('multi',plots)
        for isauto = [0,1]
            if isauto
                folder = 'multipanel-auto';
                title_suffix = ' (autocolored)';
            else
                folder = 'multipanel-standard';
                title_suffix = '';
            end
            suffix = 'mult';

            figure
            set(gcf,'PaperPosition',[0,0,9,3.5])
            tlo = tiledlayout(3,3);
            title(tlo,[plot_title,title_suffix],'FontSize',fts,'FontName',ftn,'FontWeight','bold','Interpreter','none')

            nexttile
            pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(v3_p(:,mlon_rid,:)))
            title(['N-S Ion Flow (',num2str(mlon_rac_p),'°)'])
            xlabel(mlat_label)
            ylabel(alt_label)
            colormap(gca,colorcet(clm.v))
            clb = colorbar;
            clb.Label.String = ['v_N [',units.v,']'];
            clb.Ruler.TickLabelFormat = clb_fmt;
            clb.Ruler.Exponent = clb_exp;
            if isauto
                clim(v3_range_p)
            else
                clim(v_range_hard_p)
            end

            nexttile(2,[2,1])
            pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(log10(ne_p(alt_rid,:,:))))
            title(['Density (',num2str(alt_rac_p),' ',units.x,')'])
            xlabel(mlon_label)
            ylabel(mlat_label)
            colormap(gca,colorcet(clm.n))
            clb = colorbar;
            clb.Label.String = ['log_{10} n_e [',units.n,']'];
            clb.Ruler.TickLabelFormat = clb_fmt;
            clb.Ruler.Exponent = clb_exp;
            if isauto
                clim(log10(n_range_p))
            else
                clim(log10(n_range_hard_p))
            end
            xline(mlon_rac_p,'r--')

            nexttile
            pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(jP_p(:,mlon_rid,:)))
            title(['Pedersen Current (',num2str(mlon_rac_p),'°)'])
            xlabel(mlat_label)
            ylabel(alt_label)
            colormap(gca,colorcet(clm.j))
            clb = colorbar;
            clb.Label.String = ['j_P [',units.j,']'];
            clb.Ruler.TickLabelFormat = clb_fmt;
            clb.Ruler.Exponent = clb_exp;
            if isauto
                clim(jP_range_p)
            else
                clim(j_range_hard_p)
            end

            nexttile
            pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(v2_p(:,mlon_rid,:)))
            title(['E-W Ion Flow (',num2str(mlon_rac_p),'°)'])
            xlabel(mlat_label)
            ylabel(alt_label)
            colormap(gca,colorcet(clm.v))
            clb = colorbar;
            clb.Label.String = ['v_E [',units.v,']'];
            clb.Ruler.TickLabelFormat = clb_fmt;
            clb.Ruler.Exponent = clb_exp;
            if isauto
                clim(v2_range_p)
            else
                clim(v_range_hard_p)
            end

            nexttile
            pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(jH_p(:,mlon_rid,:)))
            title(['Hall Current (',num2str(mlon_rac_p),'°)'])
            xlabel(mlat_label)
            ylabel(alt_label)
            colormap(gca,colorcet(clm.j))
            clb = colorbar;
            clb.Label.String = ['j_H [',units.j,']'];
            clb.Ruler.TickLabelFormat = clb_fmt;
            clb.Ruler.Exponent = clb_exp;
            if isauto
                clim(jH_range_p)
            else
                clim(j_range_hard_p)
            end

            nexttile
            pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),hsv_mlon)
            title(['Ion Flow (',num2str(mlon_rac_p),'°)'])
            xlabel(mlat_label)
            ylabel(alt_label)
            colormap(gca,hsv_mlon_map)
            clb = colorbar;
            colormap(clb,hsv_map_clb)
            clb.Limits = [0,1];
            clb.Ticks = [0,1/4,1/2,3/4,1];
            clb.TickLabels = {'W','S','E','N','W'};
            clb.Label.String = ['Sat. at ',num2str(hsv_sat_p),' ',units.v];
            clim([0,1])
            ylim([min(ALT_p(:)),alt_hsv_p])

            nexttile
            pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(log10(ne_p(:,mlon_rid,:))))
            title(['Density (',num2str(mlon_rac_p),'°)'])
            xlabel(mlat_label)
            ylabel(alt_label)
            colormap(gca,colorcet(clm.n))
            clb = colorbar;
            clb.Label.String = ['log_{10} n_e [',units.n,']'];
            clb.Ruler.TickLabelFormat = clb_fmt;
            clb.Ruler.Exponent = clb_exp;
            if isauto
                clim(log10(n_range_p))
            else
                clim(log10(n_range_hard_p))
            end
            yline(alt_rac_p,'r--')

            nexttile
            pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(-j1_p(:,mlon_rid,:)))
            title(['FAC (',num2str(mlon_rac_p),'°)'])
            xlabel(mlat_label)
            ylabel(alt_label)
            colormap(gca,colorcet(clm.j))
            clb = colorbar;
            clb.Label.String = ['j_{||} [',units.j,']'];
            clb.Ruler.TickLabelFormat = clb_fmt;
            clb.Ruler.Exponent = clb_exp;
            if isauto
                clim(j1_range_p)
            else
                clim(j_range_hard_p)
            end

            if ~exist(fullfile(direc,'plots',folder),'dir')
                mkdir(direc,fullfile('plots',folder));
            end
            filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
            disp(['Saving: ',filename])
            saveas(gcf,filename)
            close all
        end
    end

    % temperature plot
    if ismember('temp',plots)
        folder = 'temperature';
        suffix = 'temp';

        figure
        set(gcf,'PaperPosition',[0,0,6.5,4.5])
        tlo = tiledlayout(3,2);
        title(tlo,plot_title,'FontSize',fts,'FontName',ftn,'FontWeight','bold','Interpreter','none')

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(Te_p(:,mlon_rid,:)))
        title(['Electron Temp. (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.t))
        clb = colorbar;
        clb.Label.String = ['T_e [',units.t,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(Te_range_p)
        yline(alt_rac_p,'r--')

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(Te_p(alt_rid,:,:)))
        title(['Electron Temp. (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.t))
        clb = colorbar;
        clb.Label.String = ['T_e [',units.t,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(Te_range_p)
        xline(mlon_rac_p,'r--')

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(Ti_p(:,mlon_rid,:)))
        title(['Ion Temp. (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.t))
        clb = colorbar;
        clb.Label.String = ['T_i [',units.t,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(Ti_range_p)
        yline(alt_rac_p,'r--')

        nexttile
        pcolor(squeeze(MLON(alt_rid,:,:)),squeeze(MLAT(alt_rid,:,:)),squeeze(Ti_p(alt_rid,:,:)))
        title(['Ion Temp. (',num2str(alt_rac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.t))
        clb = colorbar;
        clb.Label.String = ['T_i [',units.t,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(Ti_range_p)
        xline(mlon_rac_p,'r--')

        nexttile
        pcolor(squeeze(MLAT(:,mlon_rid,:)),squeeze(ALT_p(:,mlon_rid,:)),squeeze(joule_p(:,mlon_rid,:)))
        title(['Joule Heating (',num2str(mlon_rac_p),'°)'])
        xlabel(mlat_label)
        ylabel(alt_label)
        colormap(gca,colorcet(clm.u))
        clb = colorbar;
        clb.Label.String = ['J\cdot E [',units.u,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(u_range_p)
        yline(alt_cac_p,'r--')

        nexttile
        pcolor(squeeze(MLON(alt_cid,:,:)),squeeze(MLAT(alt_cid,:,:)),squeeze(joule_p(alt_cid,:,:)))
        title(['Joule Heating (',num2str(alt_cac_p),' ',units.x,')'])
        xlabel(mlon_label)
        ylabel(mlat_label)
        colormap(gca,colorcet(clm.u))
        clb = colorbar;
        clb.Label.String = ['J\cdot E [',units.u,']'];
        clb.Ruler.TickLabelFormat = clb_fmt;
        clb.Ruler.Exponent = clb_exp;
        clim(u_range_p)
        xline(mlon_rac_p,'r--')

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

