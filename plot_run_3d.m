function plot_run_3d(direc,plots,alt_ref,options)
arguments
    direc (1,:) char {mustBeFolder}
    plots (1,:) string
    alt_ref (1,1) double {mustBePositive}
    options.start (1,1) double {mustBeNumeric} = -1
    options.cad (1,1) double {mustBeNumeric} = -1
    options.stop (1,1) double {mustBeNumeric} = -1
end
% assertions
plot_options = ["all","fluxtubes"];
for p = plots
    assert(ismember(p,plot_options),"'"+ p + "' is not one of the plotting options.")
end
if ismember('all',plots)
    plots = plot_options;
end
if direc(end)==filesep
    direc = direc(1:end-1);
end

%% hard coded parameters
f_scl = 1e-3; units.f = 'kA';
j_scl = 1e+6; units.j = 'uA/m^2'; clm.j = 'D1A';
n_scl = 1e+0; units.n = 'm^{-3}';
p_scl = 1e-6; units.p = 'MW';
s_scl = 1e+0; units.s = 'S/m^3';
% u_scl = 1e+6; units.u = 'uW/m^3'; clm.u = 'L19';
x_scl = 1e-3; units.x = 'km';

xlims = [-1,1]*1500e3;
ylims = [-1,1]*400e3;
zlims = [0,1]*alt_ref*1.05;
sized = 8;
qnt = 0.99; % quantile value used to set data ranges

%% loading grid data
xg = gemini3d.read.grid(direc);
x = double(xg.x2(3:end-2));
y = double(xg.x3(3:end-2));
z = double(xg.x1(3:end-2));
[~,ubz] = min(abs(z-alt_ref));
ubz = ubz + 1; % add buffer cell
z = z(1:ubz);
lz = length(z);
dx = double(xg.dx2h);
dy = double(xg.dx3h);
dz = double(xg.dx1h(1:ubz));
[X,Y,~] = ndgrid(x,y,z);
[Xm,Ym,Zm] = meshgrid(x,y,z);
[dX,dY,dZ] = ndgrid(dx,dy,dz);
dV = dX.*dY.*dZ;

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
    %% loading simulation data
    time = datetime(ymd) + seconds(UTsec);
    time.Format = 'yyyyMMdd''T''HHmmss.SSS';
    dat = gemini3d.read.frame(direc,'time',time);
    title_time = char(dat.time);
    filename_prefix = [char(time),'UT'];
    [~,runname] = fileparts(direc);

    %% formatting simulation data
    phi = dat.Phitop;
    [Ez,Ex,Ey] = gemscr.postprocess.pot2field(xg,phi);
    [sigP,sigH,~,~] = load_conductances(direc,time,dat,cfg,xg);

    Ex = permute(Ex(1:ubz,:,:),[2,3,1]);
    Ey = permute(Ey(1:ubz,:,:),[2,3,1]);
    Ez = permute(Ez(1:ubz,:,:),[2,3,1]);
    jx = permute(dat.J2(1:ubz,:,:),[2,3,1]);
    jy = permute(dat.J3(1:ubz,:,:),[2,3,1]);
    jz = permute(dat.J1(1:ubz,:,:),[2,3,1]);
    ne = permute(dat.ne(1:ubz,:,:),[2,3,1]);
    sigP = permute(sigP(1:ubz,:,:),[2,3,1]);
    sigH = -permute(sigH(1:ubz,:,:),[2,3,1]);

    % implicit simulation data
    joule = jx.*Ex + jy.*Ey + jz.*Ez;

    % set data ranges
    j1_range = [-1,1]*quantile(abs(jz(:)),qnt);
%     u_range = [quantile(joule(:),1-qnt),quantile(joule(:),qnt)];

    %% plotting routines
    set(0,'defaultTextFontSize',0.8*sized)

%     jz_p = normalize(squeeze(jz(:,:,end)),'range');
    jz_p = -squeeze(jz(:,:,end))*j_scl;
    ne_p = normalize(ne,'range');
    sigP_p = normalize(sigP,'range'); sigH_p = normalize(sigH,'range');
    X_p = X*x_scl; Y_p = Y*x_scl;
    Xm_p = Xm*x_scl; Ym_p = Ym*x_scl; Zm_p = Zm*x_scl;

    xlims_p = xlims*x_scl;
    ylims_p = ylims*x_scl;
    zlims_p = zlims*x_scl;
    alt_ref_p = alt_ref*x_scl;

%     j1_range_p = j1_range*j_scl;

    pft0 = [200e3,0e3]*x_scl;
    pft1 = -[600e3,0e3]*x_scl;

    ne_l =   [ones(lz,1)*xlims_p(1), normalize(squeeze(log10(ne_p(end/2,end/2,:))),'range',[ylims_p(1),0]), z*x_scl];
    sigP_l = [ones(lz,1)*xlims_p(1), normalize(squeeze(log10(sigP_p(end/2,end/2,:))),'range',[ylims_p(1),0]), z*x_scl];
    sigH_l = [ones(lz,1)*xlims_p(1), normalize(squeeze(log10(sigH_p(end/2,end/2,:))),'range',[ylims_p(1),0]), z*x_scl];

    % flux tube plot
    if ismember('fluxtubes',plots)
        folder = 'plots_fluxtubes';
        suffix = 'flux';
        folder_suffix = {'iso','side','top'};
        p0 = [[800,30,300];[800,80,300];[800,130,300]]*1e3*x_scl;
        r0 = [1,1,1]*200e3*x_scl;
        r1 = [1,1,1]*20e3*x_scl;
        colors = [[0,0.8,0];[1,0.5,0];[0,0.5,1]];
        views = [[30,45];[90,0];[0,90]];
        dviews = [[10*2*(UTsec0-UTsec+tdur/2)/tdur,0];[0,0];[0,0]];

        % call fluxtube prior for speed
        tubes = {0,0,0};
        for n = 1:length(r0)
            tubes{n} = fluxtube(xg,dat,alt_ref*x_scl,p0(n,:),r0(n),r1(n),reverse=1,calculate_hull=1);
        end
        for v = 1:length(views)
            vv = views(v,:)+dviews(v,:);
            for n = 1:length(r0)
                color = colors(n,:);

                tube = tubes{n};
                verts = tube.vertices;
                c0 = tube.caps.start;
                c1 = tube.caps.end;
                flux0 = tube.flux.in*j_scl*f_scl;
                flux1 = tube.flux.out*j_scl*f_scl;
                in0 = tube.flux.area.in;
                in1 = tube.flux.area.out;
                hull = tube.hull;
                joule_tot = sum(joule.*dV.*hull,'all');

                shadow = nan(size(in0));
                shadow(in0) = 0;
                shadow(in1) = 0;

                figure(1)
                set(gcf,'PaperUnits','inches','PaperPosition',[0,0,9,6])
                title([runname,' at ',title_time,' UT'],'FontSize',sized*2,'FontWeight','bold','Interpreter','none')
                
                hold on
                pl0 = plot3(c0(:,1),c0(:,2),c0(:,3));
                pl1 = plot3(c1(:,1),c1(:,2),c1(:,3));
                pl2 = plot3(ne_l(:,1),ne_l(:,2),ne_l(:,3));
                pl3 = plot3(sigP_l(:,1),sigP_l(:,2),sigP_l(:,3));
                pl4 = plot3(sigH_l(:,1),sigH_l(:,2),sigH_l(:,3));
                stl = streamline(verts);
                pcolor(squeeze(X_p(:,:,end)),squeeze(Y_p(:,:,end)),jz_p)
%                 slice(Xm_p,Ym_p,Zm_p,permute(sigP_p,[2,1,3]),-1e3,[],[]);
                shd = pcolor(squeeze(X_p(:,:,end)),squeeze(Y_p(:,:,end)),shadow);
                hold off
                shading flat
                colormap(gca,colorcet(clm.j))
                clb = colorbar;
                clb.Label.String = ['j_{||} [',units.j,']'];
                clb.FontSize = sized;
                set(pl0,'Color',color,'LineWidth',1)
                set(pl1,'Color',color,'LineWidth',1)
                set(pl2,'Color','k','LineWidth',1)
                set(pl3,'Color','r','LineWidth',1)
                set(pl4,'Color','b','LineWidth',1)
                set(shd,'FaceAlpha',0.5)
                set(stl,'Color',[color,0.5],'LineWidth',0.5)
                set(gca,'FontSize',sized)
                view(vv)
                xlim(xlims_p)
                ylim(ylims_p)
                zlim(zlims_p)
                clim([0,1])
                pbaspect([max(xlims_p),2*max(ylims_p),2*max(zlims_p)])
                xlabel(['East [',units.x,']'],'FontSize',sized)
                ylabel(['North [',units.x,']'],'FontSize',sized)
                zlabel(['Up [',units.x,']'],'FontSize',sized)
                text(c0(1,1)+pft0(1),c0(1,2)+pft0(2),1.04*alt_ref_p,[...
                    num2str(flux0,3),' ',units.f...
                    ,' ',num2str(joule_tot*p_scl,3),' ',units.p...
                    ],'Color',color)
                text(c1(1,1)+pft1(1),c1(1,2)+pft1(2),1.04*alt_ref_p,[num2str(flux1,3),' ',units.f],'Color',color)
            end
            full_folder = [folder,'_',cell2mat(folder_suffix(v))];
            if ~exist(fullfile(direc,full_folder),'dir')
                mkdir(direc,full_folder);
            end
            saveas(gcf,fullfile(direc,full_folder,[filename_prefix,'_',suffix,'.png']))
            close all
        end
    end
end
end