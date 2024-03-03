% load simulation structures
direc = '\\dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\isinglass_80\';
if not(ispc)
    direc = strrep(direc,'\','/');
end

cfg = gemini3d.read.config(direc);
if not(exist('xg','var'))
    xg = gemini3d.read.grid(direc);
    xg = jules.tools.shrink(xg);
end
if not(exist('dat','var'))
    dat = gemini3d.read.frame(direc,'time',cfg.times(end) ...
        ,'vars',["J1","J2","J3","ne","v2","v3"]);
end
fprintf('Simulation data loaded.\n')

% plotting parameters
save_plot = true;
sffx = '';
fntn = 'Arial'; % font name
fnts = 20; % font size
linw = 1.5; % line width
angl = [-25, 30]; % view angle (°)
qntl = 0.95; % colorbar range quantile
pprw = 13; % paper width (inches)
pprh = 9; % paper height (inches)
clbh = 0.43; % colorbar height (relative)
clc0 = [0, 0, 0]; % start curve color (rgb)
clc1 = [0, 0, 1]; % end curve color (rgb)
clef = [1, 0, 1]; % electric field color (rgb)
offs = 0.5; % projection line offset (km)
zoom = 1.1; % camera zoom
panx = 1.0; % pan right (°)
pany = 0.0; % pan up (°)

scl.x = 1e-3; unt.x = 'km';
scl.e = 1e3; unt.e = 'mV/m';
scl.f = 1e-3; unt.f = 'kA';
scl.n = 1e0; unt.n = 'm^{-3}'; clm.n = 'L9';
scl.j = 1e6; unt.j = 'uA/m^2'; clm.j = 'D1A';
scl.v = 1e3; unt.v = 'km/s';

lim.x = [-1, 1]*105;
lim.y = [-57, 15];
lim.z = [80, 205];

% unpack grid
x = double(xg.x2(3:end-2)*scl.x);
y = double(xg.x3(3:end-2)*scl.x);
z = double(xg.x1(3:end-2)*scl.x);
[~,lbx] = min(abs(x-max(min(x),lim.x(1))));
[~,lby] = min(abs(y-max(min(y),lim.y(1))));
[~,lbz] = min(abs(z-max(min(z),lim.z(1))));
[~,ubx] = min(abs(x-min(max(x),lim.x(2))));
[~,uby] = min(abs(y-min(max(y),lim.y(2))));
[~,ubz] = min(abs(z-min(max(z),lim.z(2))));

x = x(lbx:ubx);
y = y(lby:uby);
z = z(lbz:ubz);
dx = xg.dx2h(lbx:ubx)*scl.x;
dy = xg.dx3h(lby:uby)*scl.x;
dz = xg.dx1h(lbz:ubz)*scl.x;
[Xm,Ym,Zm] = meshgrid(x,y,z);
[dXm,dYm] = meshgrid(dx,dy);
Bmag = mean(xg.Bmag(:));

reverse_test_tolerance = min([dx;dy;dz]);
ar = [range(x),2*range(y),range(z)];

% tube parameters
tube_list = 1:3;
p0 = [ ...
    [60, -41, 200]; ...
    [10, y(end), 100]; ...
    [-74, -42, 200]; ...
    [-60, -20, 200];
    ];
v0 = [[3,1,0]; [1,0,0]; [1,0,0]; [1,0,0]];
v1 = [[0,1,0]; [0,0,1]; [0,1,0]; [0,1,0]];
r0 = [20, 20, 20, 20];
r1 = [8, 10, 4, 6];
colors = [ ...
    [1.0, 0.0, 0.0]; ...
    [1.0, 0.5, 0.0]; ...
    [0.0, 0.5, 0.0]; ...
    [1.0, 0.5, 0.7]; ...
    ];
res = [500, 400, 400, 400];
rev = [0, 1, 0, 0];
split_factor = [5, 10, 10, 10];
do_inline = [0, 0, 1, 1];
outline_axis = [3, 2, 3, 3];
outline_res = 10;

% unpack data
j1 = squeeze(dat.J1(ubz,lbx:ubx,lby:uby))'*scl.j; % A/km^2 = uA/m^2
ne = log10(squeeze(dat.ne(lbz:ubz,length(x)/2,lby:uby))')*scl.n; % m^-3
v2 = permute(squeeze(dat.v2(ubz-1:ubz,lbx:ubx,lby:uby)),[3,2,1])*scl.v; % km/s
v3 = permute(squeeze(dat.v3(ubz-1:ubz,lbx:ubx,lby:uby)),[3,2,1])*scl.v; % km/s
E2 = v3*Bmag*scl.e/scl.v; % mV/m
E3 = -v2*Bmag; % mV/m
E2(:,:,1) = nan;
E3(:,:,1) = nan;

% generate current flux tubes
for n = tube_list
    tubes.(char(64+n)) = jules.tools.current_flux_tube(xg, dat ...
        , p0(n,:), r0(n), r1(n), v0=v0(n,:), v1=v1(n,:) ...
        , res = res(n), rev = rev(n), split_factor = split_factor(n) ...
        , outline_axis = outline_axis(n), outline_res = outline_res ...
        , xlims = lim.x, ylims = lim.y, zlims = lim.z ...
        );
    fprintf('Tube %s done.\n',char(64+n))
end

%% plotting
close all
reset(0)
set(0,'defaultFigurePaperUnits','inches')
set(0,'defaultSurfaceEdgeColor','flat')
set(0,'defaultLineLineWidth',linw)
set(0,'defaultQuiverLineWidth',linw)
jules.tools.setall(0,'FontName',fntn)
jules.tools.setall(0,'FontSize',fnts)
jules.tools.setall(0,'Multiplier',1)
colorcet = @jules.tools.colorcet;

lim.j = [-1,1]*quantile(abs(j1(:)),qntl);
lim.n = [quantile(abs(ne(:)),1-qntl),quantile(abs(ne(:)),qntl)];

close all
figure
set(gcf,'Position',[50,100,1400,750])
set(gcf,'PaperUnits','inches','PaperPosition',[0,0,pprw,pprh])
axj = axes;
axn = axes;
axa = [axj,axn];
set(axj,'Color','none','XColor','none','YColor','none','ZColor','none')
set(axn,'Color','none','XColor','none','YColor','none','ZColor','none')
hold on

% current flux tubes
in = false(size(j1));
out = false(size(j1));
for n = tube_list
    color = colors(n,:);
    tube = tubes.(char(64+n));
    verts = tube.vertices;
    c0 = tube.caps.start;
    c1 = tube.caps.end;
    in = in | tube.flux.area.in;
    out = out | tube.flux.area.out;
    flux0 = tube.flux.in*scl.f;
    flux1 = tube.flux.out*scl.f;
    outline = tube.outline;
    inline = tube.inline;
    fprintf('Tube %s: influx = %.2f %s, outflux = %.2f %s, ratio = %.2f.\n' ...
        , char(64+n), flux0, unt.f, flux1, unt.f, flux0/flux1)

    stl = streamline(verts);
    plot3(c0(:,1),c0(:,2),c0(:,3),'Color',clc0);
    plot3(c0(:,1),c0(:,2),c0(:,3)*0+z(1)+offs,'Color',clc0);
    % plot3(c0(:,1)*0+x(end)-off,c0(:,2),c0(:,3),'Color',cl0);
    if ~iscell(c1)
        plot3(c1(:,1),c1(:,2),c1(:,3),'Color',clc1);
        plot3(c1(:,1),c1(:,2),c1(:,3)*0+z(1)+offs,'Color',clc1);
        % plot3(c1(:,1)*0+x(end)-off,c1(:,2),c1(:,3),'Color',cl1);
    else
        for i=1:length(c1)
            c = cell2mat(c1(i));
            plot3(c(:,1),c(:,2),c(:,3),'Color',clc1);
            plot3(c(:,1),c(:,2),c(:,3)*0+z(1),'Color',clc1);
            % plot3(c(:,1)*0+x(end)-off,c(:,2),c(:,3),'Color',cl1);
        end
    end
    plot3(outline(:,1)*0+x(end)-offs,outline(:,2),outline(:,3),'Color',color)
    if do_inline(n)
        plot3(inline(:,1)*0+x(end)-offs,inline(:,2),inline(:,3),'Color',color)
    end
    set(stl,'Color',[color,0.1],'LineWidth',linw/2)
end

% field-aligned current slice
j1(in | out) = nan;
j1_slice = repmat(j1,[1,1,length(z)]);
slice(axj,Xm,Ym,Zm,-j1_slice,[],[],z(1))
colormap(axj,colorcet(clm.j))
clim(axj,lim.j)
clb = colorbar(axj);
clb.Label.String = sprintf('j_{||} (%s)',unt.j);
clb.Position = [0.89,clbh+2*(1-2*clbh)/3,0.015,clbh];

% electron density slice
ne_slice = permute(repmat(ne,[1,1,length(x)]),[1,3,2]);
slice(axn,Xm,Ym,Zm,ne_slice,x(end),[],[])
colormap(axn,colorcet(clm.n))
clim(axn,lim.n)
clb = colorbar(axn);
clb.Label.String = sprintf('log_{10} n_e (%s)',unt.n);
clb.Position = [0.89,(1-2*clbh)/3,0.015,clbh];

% electric field
n_qx = 3;
qcy = 20;
qx_ids = round(linspace(1/2/n_qx,1-1/2/n_qx,n_qx)*length(x));
Xm_q = Xm(qcy:2*qcy:end-qcy+1,qx_ids,1:2);
Ym_q = Ym(qcy:2*qcy:end-qcy+1,qx_ids,1:2);
Zm_q = Zm(qcy:2*qcy:end-qcy+1,qx_ids,1:2);
E2_q = E2(qcy:2*qcy:end-qcy+1,qx_ids,:);
E3_q = E3(qcy:2*qcy:end-qcy+1,qx_ids,:);
E1_q = zeros(size(E2_q));
quiver3(Xm_q,Ym_q,Zm_q,E2_q,E3_q,E1_q,0.8,'-o','Color',clef,'MarkerSize',2)
quiver3(-95,y(end),185,20,0,0,0.8,'-o','Color',clef,'MarkerSize',2)
text(-95+20,y(end),185,sprintf('20 %s',unt.e),'Color',clef,'FontSize',fnts*0.8)

view(axa,angl)
xlim(axa,1.03*lim.x); ylim(axa,1.01*lim.y); zlim(axa,lim.z.*[1,1.01])
pbaspect(axj,ar); pbaspect(axn,ar)

xlabel(axj,sprintf('Mag. east (%s)',unt.x))
ylabel(axj,sprintf('Mag. north (%s)',unt.x))
zlabel(axj,sprintf('Mag. up (%s)',unt.x))

camzoom(axj,zoom); camzoom(axn,zoom)
campan(axj,panx,pany); campan(axn,panx,pany); 

% save figure
if save_plot
    if ~isempty(sffx)
        sffx = ['_',sffx];
    end
    filename = fullfile('plots','paper0',['simulation-results',sffx,'.png']);
    fprintf('Saving %s\n',filename)
    saveas(gcf,filename)
    close all
end
