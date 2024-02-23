function plot_run_3d(direc,plots,alt_ref,lon_ref,opts)
arguments
    direc (1,:) char {mustBeFolder}
    plots (1,:) string
    alt_ref (1,1) double {mustBePositive}
    lon_ref (1,1) double {mustBeNumeric}
    opts.start (1,1) double {mustBeNumeric} = -1
    opts.cad (1,1) double {mustBeNumeric} = -1
    opts.stop (1,1) double {mustBeNumeric} = -1
    opts.cfg (:,:) struct = struct
    opts.xg (:,:) struct = struct
    opts.dat (:,:) struct = struct
    opts.plot_type (1,1) int8 = 1
    opts.zoom (1,1) double = 1
end

doAnnotate = false;
paper_w = [13,5,8];
paper_h = [9,5,5];

% assertions
plot_options = ["all","fluxtubes"];
for p = plots
    assert(ismember(p,plot_options),"'"+ p + "' is not one of the plotting options.")
end
if ismember('all',plots)
    plots = plot_options;
end
if any(direc(end)=='/\')
    direc = direc(1:end-1);
end

%% hard coded parameters
scl.f = 1e-3; units.f = 'kA';
scl.j = 1e+6; units.j = 'uA/m^2';   clm.j = 'D1A';
scl.n = 1e+0; units.n = 'm^{-3}';   clm.n = 'L9';
scl.p = 1e-6; units.p = 'MW';
scl.s = 1e+0; units.s = 'S/m';      clm.s = 'L18';
scl.x = 1e-3; units.x = 'km';

colorcet = @jules.tools.colorcet;

zmin = 80e3;
% xlims = [-1,1]*1050e3;
% xlims = [-1,1]*1300e3;
xlims = [-1,1]*103e3;
ylims = [-57,0]*1e3;
% ylims = [-190,360]*1e3;
% ylims2 = [-120,290]*1e3;
% ylims2 = ylims;
zlims = [zmin,alt_ref*1.05];
% qnt = 0.99; % quantile value used to set data ranges
qnt = 0.95;
fts = 20; %17 % fontsize
ftn = 'Arial';
if ispc
    lnw = 1;
else
    lnw = 0.5;
end
% clb_fmt = '%+ 5.1f'; % colorbar ticklabel format
% clb_fmt = '%3.0f';
clb_exp = 0; % force no colorbar exponents

%% loading grid data
if isempty(fields(opts.xg))
    xg = gemini3d.read.grid(direc);
else
    xg = opts.xg;
end
x = double(xg.x2(3:end-2));
y = double(xg.x3(3:end-2));
z = double(xg.x1(3:end-2));
[~,lbx] = min(abs(x-xlims(1))); [~,ubx] = min(abs(x-xlims(2)));
[~,lby] = min(abs(y-ylims(1))); [~,uby] = min(abs(y-ylims(2)));
[~,ubz] = min(abs(z-alt_ref));
ubz = ubz + 1; % add buffer cell
x = x(lbx:ubx); y = y(lby:uby); z = z(1:ubz);
lz = length(z);
% dx = double(xg.dx2h);
% dy = double(xg.dx3h);
% dz = double(xg.dx1h(1:ubz));
[Xm,Ym,Zm] = meshgrid(x,y,z);
% [dX,dY,dZ] = ndgrid(dx,dy,dz);
% dV = dX.*dY.*dZ;

%% loading configuration data
if isempty(fields(opts.cfg))
    cfg = gemini3d.read.config(direc);
else
    cfg = opts.cfg;
end
ymd = cfg.ymd;
UTsec0 = cfg.UTsec0;
tdur = cfg.tdur;
dtout = cfg.dtout;

%% setting time boundaries
if opts.start < 0
    start = 0;
else
    start = opts.start;
end
if opts.cad <= 0
    cad = dtout;
else
    cad = opts.cad;
end
if opts.stop < start
    stop = tdur;
else
    stop = opts.stop;
end

%% main loop
for UTsec = UTsec0+start:cad:UTsec0+stop
    %% loading simulation data
    time = datetime(ymd) + seconds(UTsec);
    time.Format = 'yyyyMMdd''T''HHmmss.SSS';
    if isempty(fields(opts.dat))
        dat = gemini3d.read.frame(direc,'time',time);
    else
        dat = opts.dat;
    end
    title_time = char(dat.time);
    % filename_prefix = [char(time),'UT'];
    filename_prefix = char(gemini3d.datelab(time));
    [~,runname] = fileparts(direc);

    %% formatting simulation data
    %     phi = dat.Phitop;
    %     [Ez,Ex,Ey] = gemscr.postprocess.pot2field(xg,phi);

    %     Ex = permute(Ex(1:ubz,:,:),[2,3,1]);
    %     Ey = permute(Ey(1:ubz,:,:),[2,3,1]);
    %     Ez = permute(Ez(1:ubz,:,:),[2,3,1]);
    %     jx = permute(dat.J2(1:ubz,:,:),[2,3,1]);
    %     jy = permute(dat.J3(1:ubz,:,:),[2,3,1]);
    jz = permute(dat.J1(1:ubz,lbx:ubx,lby:uby),[2,3,1]);
    ne = permute(dat.ne(1:ubz,lbx:ubx,lby:uby),[2,3,1]);

    % implicit simulation data
    %     joule = jx.*Ex + jy.*Ey + jz.*Ez;

    % set data ranges
    j1_range = [-1,1]*quantile(abs(jz(:,:,end)),qnt,'all');

    %% plotting routines
    set(0,'defaultTextFontSize',fts)

    jz_p = -jz*scl.j;
    ne_p = ne*scl.n;
    Xm_p = Xm*scl.x; Ym_p = Ym*scl.x; Zm_p = Zm*scl.x;

    xlims_p = xlims*scl.x;
    ylims_p = ylims*scl.x;
    zlims_p = zlims*scl.x;
    lon_ref_p = lon_ref*scl.x;
    alt_ref_p = alt_ref*scl.x;

    j1_range_p = j1_range*scl.j;

    % flux tube plot
    if ismember('fluxtubes',plots)
        folder = ['fluxtubes',num2str(opts.plot_type)];
        suffix = 'flux';
        folder_suffix = {'iso','side','top'};

        if opts.plot_type == 1
            ntubes = 3;
            p0 = [[lon_ref_p,40,alt_ref_p];[lon_ref_p,110,alt_ref_p];[lon_ref_p,180,alt_ref_p]]*1e3*scl.x;
            v0 = repmat([1,0,0],ntubes,1);
            v1 = repmat([0,1,0],ntubes,1);
            r0 = ones(1,ntubes)*200e3*scl.x;
            r1 = ones(1,ntubes)*20e3*scl.x;
            colors = [...
                [1.0, 0.5, 0.0];...
                [0.2, 0.2, 0.8];...
                [0.0, 0.5, 0.0];...
                ];
            res = [64,64,64];
            rev = ones(1,ntubes);
        elseif opts.plot_type == 2
            ntubes = 3;
            p0 = [[lon_ref_p,40,300];[lon_ref_p,160,300];[1000,100,115]]*1e3*scl.x;
            v0 = [[1,0,0];[1,0,0];[0,0,1]];
            v1 = repmat([0,1,0],ntubes,1);
            r0 = [200,200,10]*1e3*scl.x;
            r1 = [20,20,45]*1e3*scl.x;
            colors = [...
                [1.0, 0.5, 0.0];...
                [0.2, 0.2, 0.8];...
                [0.0, 0.0, 0.0]...
                ];
            res = [64,64,32];
            rev = ones(1,ntubes);
        elseif opts.plot_type == 3
            ntubes = 10;
            p0 = [-400+30*(0:10:40); ones(1,ntubes/2)*7; 125+(0:10:40)]'*1e3*scl.x;
            p0 = repmat(p0,2,1);
            v0 = repmat([1,0,0],ntubes,1);
            v1 = repmat([0,0,1],ntubes,1);
            r0 = ones(1,ntubes)*100e3*scl.x;
            r1 = ones(1,ntubes)*5e3*scl.x;
            %             colors = [(1-(0:10:40)/40); (0:10:40)/40; zeros(1,5)]';
            colors = [[1,0,0];[1,0,1];[0,0,1];[0,1,1];[0,1,0]];
            colors = repmat(colors,2,1);
            res = ones(1,ntubes)*32;
            rev = [ones(1,ntubes/2),zeros(1,ntubes/2)];
        elseif opts.plot_type == 4
            ntubes = 3;
            p0 = [[lon_ref_p,5,200];[lon_ref_p,10,200];[lon_ref_p,15,200]]*1e3*scl.x;
            v0 = repmat([1,0,0],ntubes,1);
            v1 = repmat([0,1,0],ntubes,1);
            r0 = ones(1,ntubes)*2e3*scl.x;
            r1 = ones(1,ntubes)*2e3*scl.x;
            colors = [...
                [1.0, 0.5, 0.0];...
                [0.2, 0.2, 0.8];...
                [0.0, 0.5, 0.0];...
                ];
            res = [64,64,64];
            rev = ones(1,ntubes);
        elseif opts.plot_type == 5
            ntubes = 3;
            p0 = [[lon_ref_p,40+50,300];[lon_ref_p,110+50,300];[lon_ref_p,180+50,300]]*1e3*scl.x;
            v0 = repmat([1,0,0],ntubes,1);
            v1 = repmat([0,1,0],ntubes,1);
            r0 = ones(1,ntubes)*200e3*scl.x;
            r1 = ones(1,ntubes)*20e3*scl.x;
            colors = [...
                [1.0, 0.5, 0.0];...
                [0.2, 0.2, 0.8];...
                [0.0, 0.5, 0.0];...
                ];
            res = [64,64,64];
            rev = ones(1,ntubes);
        elseif opts.plot_type == 6
            ntubes = 3;
            p0 = [[lon_ref_p,40,alt_ref/1e3];[lon_ref_p,110,alt_ref/1e3];[lon_ref_p,180,alt_ref/1e3]]*1e3*scl.x;
            v0 = repmat([1,0,0],ntubes,1);
            v1 = repmat([0,1,0],ntubes,1);
            r0 = ones(1,ntubes)*200e3*scl.x;
            r1 = ones(1,ntubes)*20e3*scl.x;
            colors = [...
                [1.0, 0.5, 0.0];...
                [0.2, 0.2, 0.8];...
                [0.0, 0.5, 0.0];...
                ];
            res = [64,64,64];
            rev = ones(1,ntubes);
        elseif opts.plot_type == 7
            ntubes = 3;
            p0 = [[lon_ref_p,40,alt_ref/1e3];[lon_ref_p,100,alt_ref/1e3];[lon_ref_p,195,alt_ref/1e3]]*1e3*scl.x;
            v0 = repmat([1,0,0],ntubes,1);
            v1 = repmat([0,1,0],ntubes,1);
            r0 = ones(1,ntubes)*200e3*scl.x;
            r1 = ones(1,ntubes)*20e3*scl.x;
            colors = [...
                [1.0, 0.5, 0.0];...
                [0.2, 0.2, 0.8];...
                [0.0, 0.5, 0.0];...
                ];
            res = [64,64,64];
            rev = ones(1,ntubes);
        elseif opts.plot_type == 8
            ntubes = 4;
            % p0 = [[57,-23,alt_ref/1e3];[15,-43,alt_ref/1e3];[-55,-31,alt_ref/1e3]]*1e3*scl.x;
            p0 = [ ...
                [51,-23,alt_ref/1e3]; ...
                [22,-38,alt_ref/1e3]; ...
                [22,-45.5,alt_ref/1e3]; ...
                [-60,-31,alt_ref/1e3] ...
                ]*1e3*scl.x;
            v0 = [[4,1,0];[10,1,0];[10,1,0];[-13,1,0]];
            v1 = repmat([0,1,0],ntubes,1);
            r0 = [20,10,18,20]*1e3*scl.x;
            r1 = [3,2,3,3]*1e3*scl.x;
            colors = [...
                [1.0, 0.5, 0.0];...
                [0.0, 0.5, 0.0];...
                [0.0, 0.5, 0.0];...
                [1.0, 0.0, 0.0];...
                ];
            res = [200,100,100,200]*2;
            rev = [1,0,0,1];
        end

%         views = [[30,45];[90,0];[0,90]];
        views = [[225-10-180-8,38];[90,0];[0,90]];
%         views = [[-110,-30];[90,0];[0,90]];
        dviews = [[10*2*(UTsec0-UTsec+tdur/2)/tdur,0];[0,0];[0,0]];
%         dviews = [[10*2*(UTsec0-36000+150/2)/150,0];[0,0];[0,0]];
%         paper_w = [8,11,9.5];
%         paper_w = [14,5.3,5.7];
        fluxes = zeros(2,ntubes);
        joule_heatings = zeros(1,ntubes);

        % call fluxtube prior for speed
        tubes = struct;
        for n = 1:ntubes
            tubes.(char(64+n)) = jules.tools.fluxtube(xg,dat,alt_ref*scl.x...
                ,p0(n,:),r0(n),r1(n),v0=v0(n,:),v1=v1(n,:)...
                ,reverse=rev(n),calculate_hull=0,res=res(n)...
                ,xlims = xlims_p,ylims = ylims_p);
        end
        for v = 1%:length(views)
            vv = views(v,:)+dviews(v,:);

            figure(v)
            set(gcf,'PaperUnits','inches','PaperPosition',[0,0,paper_w(v),paper_h(v)])
            title([runname,' at ',title_time,' UT'],'FontSize',fts*2,'FontWeight','bold','Interpreter','none')
            t = tiledlayout(1,1,'TileSpacing','compact');
            axj = axes(t); %#ok<LAXES>
            axn = axes(t); %#ok<LAXES>
            axt = axes(t); %#ok<LAXES>
            axa = [axj,axn,axt];

            set(axa,'FontSize',fts)
            set(axa,'FontName',ftn)
            set(axj,'Color','none','XColor','none','YColor','none','ZColor','none','ZTick',alt_ref_p)
            set(axn,'Color','none','XColor','none','YColor','none','ZColor','none','XTick',lon_ref_p)
            set(axt,'Color','none','XGrid','on','YGrid','on','ZGrid','on')

            xlim(axj,xlims_p)
            xlim(axn,xlims_p - (lon_ref_p+xlims_p(1)))
            xlim(axt,xlims_p)
            
            ylim(axj,ylims_p)
            ylim(axn,ylims_p)
            ylim(axt,ylims_p)

            zlim(axj,zlims_p + (alt_ref-zmin)*scl.x)
            zlim(axn,zlims_p)
            zlim(axt,zlims_p)
            
            view(axa,vv)
%             ar = [0.4*range(xlims_p)*0.7,range(ylims_p),2*range(zlims_p)];
            ar = [range(xlims_p),range(ylims_p)*2,range(zlims_p)];
            pbaspect(axj,ar)
            pbaspect(axn,ar)
            pbaspect(axt,ar)
            hold(axa,'on')

            slice(axj,Xm_p,Ym_p,Zm_p,permute(jz_p,[2,1,3]),[],[],alt_ref_p);
            colormap(axj,colorcet(clm.j))
            shading(axj,'flat')
            clim(axj,j1_range_p)
            if v == 1%2
                clb = colorbar(axj);
                clb.Label.String = ['j_{||} (',units.j,')'];
                clb.FontSize = fts;
                clb.Position = [0.89,0.12,0.015,0.41];
                % clb.Ruler.TickLabelFormat = clb_fmt;
                % clb.Ruler.Exponent = clb_exp;
            end
            
            slice(axn,Xm_p,Ym_p,Zm_p,permute(log10(ne_p),[2,1,3]),lon_ref_p,[],[]);
            colormap(axn,colorcet(clm.n))
            clim(axn,[9.8,11.9])
            shading(axn,'flat')
            if v == 1%2
                clb = colorbar(axn);
                clb.Label.String = ['log n_e (',units.n,')'];
                clb.FontSize = fts;
                clb.Position = [0.89,0.57,0.015,0.41];
%                 clb.Ruler.TickLabelFormat = clb_fmt;
                clb.Ruler.Exponent = clb_exp;
            end

            for n = 1:ntubes
                color = colors(n,:);
                tube = tubes.(char(64+n));
                verts = tube.vertices;
                c0 = tube.caps.start;
                c1 = tube.caps.end;
                fluxes(1,n) = tube.flux.in*scl.j*scl.f;
                fluxes(2,n) = tube.flux.out*scl.j*scl.f;
                in0 = tube.flux.area.in;
                in1 = tube.flux.area.out;
                % hull = tube.hull;
                % joule_heatings(n) = sum(joule.*dV.*hull,'all');

                shadow = nan(size(in0));
                shadow(in0) = 0;
                shadow(in1) = 0;
                shadow = repmat(shadow,[1,1,lz]);

                pl0 = plot3(axt,c0(:,1),c0(:,2),c0(:,3));
                pl1 = plot3(axt,c1(:,1),c1(:,2),c1(:,3));
                if all(range(c0(:,3))<1)
                    pl0s = plot3(axt,c0(:,1),c0(:,2),ones(size(c0(:,3)))*zmin*scl.x);
                end
                if all(range(c1(:,3))<1)
                    pl1s = plot3(axt,c1(:,1),c1(:,2),ones(size(c1(:,3)))*zmin*scl.x);
                end
                stl = streamline(axt,verts);
                if any(not(isnan(shadow)),'all')
                    shd = slice(axj,Xm_p,Ym_p,Zm_p,permute(shadow,[2,1,3]),[],[],alt_ref_p);
                end
%                 plot3(axt,[1,1,1]*lon_ref_p,[ylims_p,ylims_p(2)],[zlims_p(1),zlims_p],'--','Color',[0,1,0]);
%                 plot3(axt,[xlims_p(1),xlims_p],[ylims_p,ylims_p(2)],[1,1,1]*alt_ref_p,'--','Color',[0,0,1]);

                shading(axj,'flat')
                set(pl0,'Color','k','LineWidth',lnw*2)
                set(pl1,'Color','b','LineWidth',lnw*2)
                set(pl0s,'Color','k','LineWidth',lnw*2,'LineStyle',':')
                set(pl1s,'Color','b','LineWidth',lnw*2,'LineStyle',':')
                set(shd,'FaceAlpha',0.5)
                set(stl,'Color',[color,0.2],'LineWidth',lnw)
            end

            xlabel(['Mag. east (',units.x,')'],'FontSize',fts)
            ylabel(['Mag. north (',units.x,')'],'FontSize',fts)
            zlabel(['Mag. up (',units.x,')'],'FontSize',fts)

            flux_strings = cell(2,ntubes);
            heat_strings = cell(1,ntubes);
            for n = 1:ntubes
                flux_strings{1,n} = ['\color[rgb]{',num2str(colors(n,:)),'}'...
                    ,num2str(fluxes(1,n),'%3.1f'),' ',units.f,'  ','\color{black}'];
                flux_strings{2,n} = ['\color[rgb]{',num2str(colors(n,:)),'}'...
                    ,num2str(fluxes(2,n),'%3.1f'),' ',units.f,'  ','\color{black}'];
                heat_strings{n} = ['\color[rgb]{',num2str(colors(n,:)),'}'...
                    ,num2str(joule_heatings(n)*scl.p,'%3.1f'),' ',units.p,'  ','\color{black}'];
            end
            if v == 1 && doAnnotate
                annotation('textbox',[0.52,0.97,0.01,0.01],'String',[...
                    'Current In:  ',cell2mat(flux_strings(1,:)),newline,...
                    'Current Out:  ',cell2mat(flux_strings(2,:)),newline,...
%                     'Joule Heating:  ',cell2mat(heat_strings),...
                    ],'FitBoxToText','on','EdgeColor','none','BackgroundColor','none','FontSize',18,'FontName',ftn) %#ok<*UNRCH>
            end
            if v == 1
                pan_x = 0.0;
                pan_y = 0.7; % >0 moves image down
                zoom = opts.zoom;
                campan(axj,pan_y,pan_x); campan(axn,pan_y,pan_x); campan(axt,pan_y,pan_x)
                camzoom(axj,zoom); camzoom(axn,zoom); camzoom(axt,zoom)
                % annotation('textarrow',[0.12,0.075],[0.12,0.04],'String','N', ...
                %     'LineWidth',6*lnw,'HeadLength',16,'HeadWidth',16,'HeadStyle','plain')
            end

            full_folder = [folder,'_',cell2mat(folder_suffix(v))];
            if ~exist(fullfile(direc,'plots_3d',full_folder),'dir')
                mkdir(direc,fullfile('plots_3d',full_folder));
            end
            saveas(gcf,fullfile(direc,'plots_3d',full_folder,[filename_prefix,'_',suffix,'.png']))
            disp(fullfile(direc,'plots_3d',full_folder,[filename_prefix,'_',suffix,'.png']))
            close all
        end
    end
end
end
