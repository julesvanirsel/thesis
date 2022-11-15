function plot_run_3d(direc,plots,alt_ref,lon_ref,options)
arguments
    direc (1,:) char {mustBeFolder}
    plots (1,:) string
    alt_ref (1,1) double {mustBePositive}
    lon_ref (1,1) double {mustBeNumeric}
    options.start (1,1) double {mustBeNumeric} = -1
    options.cad (1,1) double {mustBeNumeric} = -1
    options.stop (1,1) double {mustBeNumeric} = -1
    options.dat (:,:) struct = struct
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
j_scl = 1e+6; units.j = 'uA/m^2';   clm.j = 'D1A';
n_scl = 1e+0; units.n = 'm^{-3}';   clm.n = 'L9';
p_scl = 1e-6; units.p = 'MW';
% s_scl = 1e+0; units.s = 'S/m';      clm.s = 'L18';
x_scl = 1e-3; units.x = 'km';

zmin = 80e3;
xlims = [-1,1]*1450e3;
ylims = [-190,360]*1e3;
zlims = [zmin,alt_ref*1.05];
sized = 10;
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
% [X,Y,~] = ndgrid(x,y,z);
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
    [Ez,Ex,Ey] = gemscr.postprocess.pot2field(xg,phi);
%     [sigP,sigH,~,~] = load_conductances(direc,time,dat,cfg,xg);

    Ex = permute(Ex(1:ubz,:,:),[2,3,1]);
    Ey = permute(Ey(1:ubz,:,:),[2,3,1]);
    Ez = permute(Ez(1:ubz,:,:),[2,3,1]);
    jx = permute(dat.J2(1:ubz,:,:),[2,3,1]);
    jy = permute(dat.J3(1:ubz,:,:),[2,3,1]);
    jz = permute(dat.J1(1:ubz,:,:),[2,3,1]);
    ne = permute(dat.ne(1:ubz,:,:),[2,3,1]);
%     sigP = permute(sigP(1:ubz,:,:),[2,3,1]);
%     sigH = -permute(sigH(1:ubz,:,:),[2,3,1]);

    % implicit simulation data
    joule = jx.*Ex + jy.*Ey + jz.*Ez;

    % set data ranges
    j1_range = [-1,1]*quantile(abs(jz(:)),qnt);

    %% plotting routines
    set(0,'defaultTextFontSize',0.8*sized)

    jz_p = -jz*j_scl;
%     sigP_p = sigP*s_scl; sigH_p = sigH*s_scl;
    ne_p = ne*n_scl;
    Xm_p = Xm*x_scl; Ym_p = Ym*x_scl; Zm_p = Zm*x_scl;

    xlims_p = xlims*x_scl;
    ylims_p = ylims*x_scl;
    zlims_p = zlims*x_scl;
    lon_ref_p = lon_ref*x_scl;
    alt_ref_p = alt_ref*x_scl;

    j1_range_p = j1_range*j_scl;

    % flux tube plot
    if ismember('fluxtubes',plots)
        folder = 'fluxtubes';
        suffix = 'flux';
        folder_suffix = {'iso','side','top','side2'};
        p0 = [[lon_ref_p,40,300];[lon_ref_p,100,300];[lon_ref_p,160,300]]*1e3*x_scl;
        r0 = [1,1,1]*400e3*x_scl;
        r1 = [1,1,1]*20e3*x_scl;
        colors = [[1, 0.5, 0];...
                  [0.2, 0.8, 0.2];...
                  [0, 0.2, 1]];
        views = [[30,45];[90,0];[0,90];[0,0]];
        dviews = [[10*2*(UTsec0-UTsec+tdur/2)/tdur,0];[0,0];[0,0];[0,0]];
        ntubes = length(r0);
%         text_pos = [[0.60, 0.70, 0.22, 0.70];...
%                     [0.65, 0.90, 0.22, 0.90];...
%                     [0.70, 0.70, 0.15, 0.35]];
%         text_spc = 0.025;
        fluxes = zeros(2,ntubes);
        joule_heatings = zeros(1,ntubes);

        % call fluxtube prior for speed
        tubes = {0,0,0};
        for n = 1:ntubes
            tubes{n} = fluxtube(xg,dat,alt_ref*x_scl,p0(n,:),r0(n),r1(n),reverse=1,calculate_hull=1,res=64);
        end
        for v = 1:length(views)
            vv = views(v,:)+dviews(v,:);

            figure(v)
            set(gcf,'PaperUnits','inches','PaperPosition',[0,0,12,6])
            title([runname,' at ',title_time,' UT'],'FontSize',sized*2,'FontWeight','bold','Interpreter','none')
            t = tiledlayout(1,1,'TileSpacing','compact');
            axj = axes(t); %#ok<LAXES>
            axn = axes(t); %#ok<LAXES>
            axt = axes(t); %#ok<LAXES>
            axa = [axj,axn,axt];

            set(axa,...
                'FontSize',sized...
                )
            set(axj,...
                'XColor','none'...
                ,'YColor','none'...
                ,'ZColor','blue'...
                ,'ZTick',alt_ref_p...
                )
            set(axn,...
                'Color','none'...
                ,'XColor','red'...
                ,'YColor','none'...
                ,'ZColor','none'...
                ,'XTick',lon_ref_p ...
                )
            set(axt,...
                'Color','none'...
                ,'XGrid','on'...
                ,'YGrid','on'...
                ,'ZGrid','on'...
                )

            xlim(axj,xlims_p)
            xlim(axn,xlims_p + (lon_ref_p-xlims_p(1)))
            xlim(axt,xlims_p)
            ylim(axa,ylims_p)
            zlim(axj,zlims_p + (alt_ref-zmin)*x_scl)
            zlim(axn,zlims_p)
            zlim(axt,zlims_p)
            view(axa,vv)
            ar = [0.4*range(xlims_p),range(ylims_p),2*range(zlims_p)];
            pbaspect(axj,ar)
            pbaspect(axn,ar)
            pbaspect(axt,ar)
            hold(axa,'on')

            slice(axj,Xm_p,Ym_p,Zm_p,permute(jz_p,[2,1,3]),[],[],alt_ref_p);
            colormap(axj,colorcet(clm.j))
            shading(axj,'flat')
            clb = colorbar(axj);
            clb.Label.String = ['j_{||} [',units.j,']'];
            clb.FontSize = sized;
            clb.Position = [0.91,0.04,0.015,0.44];
            clim(axj,j1_range_p)

%             slice(axs,Xm_p,Ym_p,Zm_p,permute(sigP_p,[2,1,3]),xref_p,[],[]);
            slice(axn,Xm_p,Ym_p,Zm_p,permute(ne_p,[2,1,3]),lon_ref_p,[],[]);
            colormap(axn,colorcet(clm.n))
            shading(axn,'flat')
            clb = colorbar(axn);
            clb.Label.String = ['n_e [',units.n,']'];
            clb.FontSize = sized;
            clb.Position = [0.91,0.52,0.015,0.44];

            for n = 1:ntubes
                color = colors(n,:);
                tube = tubes{n};
                verts = tube.vertices;
                c0 = tube.caps.start;
                c1 = tube.caps.end;
                fluxes(1,n) = tube.flux.in*j_scl*f_scl;
                fluxes(2,n) = tube.flux.out*j_scl*f_scl;
                in0 = tube.flux.area.in;
                in1 = tube.flux.area.out;
                hull = tube.hull;
                joule_heatings(n) = sum(joule.*dV.*hull,'all');

                shadow = nan(size(in0));
                shadow(in0) = 0;
                shadow(in1) = 0;
                shadow = repmat(shadow,[1,1,lz]);

                pl0 = plot3(axt,c0(:,1),c0(:,2),c0(:,3));
                pl1 = plot3(axt,c1(:,1),c1(:,2),c1(:,3));
                stl = streamline(axt,verts);
                shd = slice(axj,Xm_p,Ym_p,Zm_p,permute(shadow,[2,1,3]),[],[],alt_ref_p);
                %                 isosurface(Xm_p,Ym_p,Zm_p,permute(hull*mean(-jz_p(:)),[2,1,3]),mean(-jz_p(:)))

                shading(axj,'flat')
                set(pl0,'Color',color,'LineWidth',2)
                set(pl1,'Color',color,'LineWidth',1)
                set(shd,'FaceAlpha',0.5)
                set(stl,'Color',[color,0.5],'LineWidth',1)

                xlabel(['East [',units.x,']'],'FontSize',sized)
                ylabel(['North [',units.x,']'],'FontSize',sized)
                zlabel(['Up [',units.x,']'],'FontSize',sized)
                %                 text(text_pos(v,1), text_pos(v,2)-(n-1)*text_spc,...
                %                     [...
                %                     num2str(flux0,3), ' ', units.f, ' ',...
                %                     num2str(joule_tot*p_scl,3), ' ', units.p...
                %                     ], 'Color', color, 'fontweight', 'bold', 'Units', 'normalized')
                %                 text(text_pos(v,3), text_pos(v,4)-(n-1)*text_spc,...
                %                     [...
                %                     num2str(flux1,3), ' ', units.f...
                %                     ], 'Color', color, 'fontweight', 'bold', 'Units', 'normalized')
%                 annotation('textbox',[text_pos(v,1),text_pos(v,2)-(n-1)*text_spc,0.01,0.01],'string',...
%                     [...
%                     ' ', num2str(flux0,3), ' ', units.f, ' - ',...
%                     num2str(joule_tot*p_scl,3), ' ', units.p...
%                     ], 'FitBoxToText','on', 'Color', color, 'EdgeColor', 'none', 'BackgroundColor', 'w' ...
%                     )
%                 annotation('textbox',[text_pos(v,3),text_pos(v,4)-(n-1)*text_spc,0.01,0.01],'string',...
%                     [...
%                     ' ', num2str(flux1,3), ' ', units.f,...
%                     ], 'FitBoxToText','on', 'Color', color, 'EdgeColor', 'none', 'BackgroundColor', 'w' ...
%                     )
            end
            
            flux_strings = cell(2,ntubes);
            heat_strings = cell(1,ntubes);
            for n = 1:ntubes
                flux_strings{1,n} = ['\color[rgb]{',num2str(colors(n,:)),'}'...
                    ,num2str(fluxes(1,n),3),' ',units.f,'  ','\color{black}'];
                flux_strings{2,n} = ['\color[rgb]{',num2str(colors(n,:)),'}'...
                    ,num2str(fluxes(2,n),3),' ',units.f,'  ','\color{black}'];
                heat_strings{n} = ['\color[rgb]{',num2str(colors(n,:)),'}'...
                    ,num2str(joule_heatings(n)*p_scl,3),' ',units.p,'  ','\color{black}'];
            end

            annotation('textbox',[0.02,0.97,0.01,0.01],'String',[...
                'Current In:  ',cell2mat(flux_strings(1,:)),newline,...
                'Current Out:  ',cell2mat(flux_strings(2,:)),newline,...
                'Joule Heating:  ',cell2mat(heat_strings),...
                ],'FitBoxToText','on','EdgeColor','none','BackgroundColor','w')

            full_folder = [folder,'_',cell2mat(folder_suffix(v))];
            if ~exist(fullfile(direc,'plots_3d',full_folder),'dir')
                mkdir(direc,fullfile('plots_3d',full_folder));
            end
            saveas(gcf,fullfile(direc,'plots_3d',full_folder,[filename_prefix,'_',suffix,'.png']))
            close all
        end
    end
end
end