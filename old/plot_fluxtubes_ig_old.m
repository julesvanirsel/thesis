direc = '\\dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\isinglass_80\';
if not(ispc)
    direc = strrep(direc,'\','/');
end
if not(exist('cfg','var'))
    cfg = gemini3d.read.config(direc);
end
if not(exist('xg','var'))
    xg = gemini3d.read.grid(direc);
    xg_tmp.x1 = xg.x1; xg_tmp.x2 = xg.x2; xg_tmp.x3 = xg.x3;
    xg_tmp.dx1h = xg.dx1h; xg_tmp.dx2h = xg.dx2h; xg_tmp.dx3h = xg.dx3h;
end
if not(exist('dat','var'))
    dat = gemini3d.read.frame(direc,'time',cfg.times(end) ...
        ,'vars',["J1","J2","J3","ne"]);
end
fprintf('Data loaded.\n')

%%
close all
tubes = plot_fluxtubes(220e3,xg=xg,dat=dat,zoom=1.25,suffix='_test',annotate=true);

%%
function tubes = plot_fluxtubes(alt_ref,opts)
arguments
    alt_ref (1,1) double {mustBePositive}
    opts.xg (:,:) struct = struct
    opts.dat (:,:) struct = struct
    opts.zoom (1,1) double = 1
    opts.suffix (1,:) char = ''
    opts.annotate (1,1) logical = false
end
tic

paper_w = [13,5,8];
paper_h = [9,5,5];

%% hard coded parameters
scl.f = 1e-3; units.f = 'kA';
scl.j = 1e+6; units.j = 'uA/m^2';   clm.j = 'D1A';
scl.n = 1e+0; units.n = 'm^{-3}';   clm.n = 'L9';
scl.x = 1e-3; units.x = 'km';

colorcet = @jules.tools.colorcet;

zmin = 80e3;
lim.x = [-1,1]*103;
lim.y = [-57,0];
lim.z = [zmin,alt_ref*1.05]*scl.x;
qnt = 0.95;
fts = 20;
ftn = 'Arial';
lnw = 1;
clb_exp = 0; % force no colorbar exponents

%% loading grid data
xg = opts.xg;
x = double(xg.x2(3:end-2));
y = double(xg.x3(3:end-2));
z = double(xg.x1(3:end-2));
[~,lbx] = min(abs(x*scl.x-lim.x(1))); [~,ubx] = min(abs(x*scl.x-lim.x(2)));
[~,lby] = min(abs(y*scl.x-lim.y(1))); [~,uby] = min(abs(y*scl.x-lim.y(2)));
[~,ubz] = min(abs(z-alt_ref));
ubz = ubz + 1; % add buffer cell
x = x(lbx:ubx); y = y(lby:uby); z = z(1:ubz);
lz = length(z);
[Xm,Ym,Zm] = meshgrid(x,y,z);

%% loading and formatting simulation data
dat = opts.dat;
jz = permute(dat.J1(1:ubz,lbx:ubx,lby:uby),[2,3,1]);
ne = permute(dat.ne(1:ubz,lbx:ubx,lby:uby),[2,3,1]);

fprintf('first: %i seconds.\n',round(toc))

%% plotting routines
jz_p = -jz*scl.j;
ne_p = ne*scl.n;
Xm_p = Xm*scl.x; Ym_p = Ym*scl.x; Zm_p = Zm*scl.x;

alt_ref_p = alt_ref*scl.x;

lim.j = [-1,1]*quantile(abs(jz(:,:,end)),qnt,'all')*scl.j;

% flux tube plot
ntubes = 3;
% p0 = [ ...
%     [51,-23,alt_ref/1e3]; ...
%     [22,-38,alt_ref/1e3]; ...
%     [22,-45.5,alt_ref/1e3]; ...
%     [-60,-31,alt_ref/1e3] ...
%     ]*1e3*scl.x;
% p0 = [ ...
%     [-50, -42, alt_ref/1e3]; ...
%     [ 20, -45, alt_ref/1e3]; ...
%     [ 80, -32, alt_ref/1e3]; ...
%     [0, -47, alt_ref/1e3];%[-90, -30, 140]; ...
%     ]*1e3*scl.x;
% v0 = [[-13,1,0]; [10,1,0]; [4,1,0]; [-13,1,0]];%[0,0,1]
% v1 = [[0,1,0]; [0,1,0]; [0,1,0]; [0,1,0];];
% r0 = [20, 20, 20, 40]*1e3*scl.x;
% r1 = [3, 3, 3, 8]*1e3*scl.x;
p0 = [ ...
    [-30, -21, alt_ref/1e3]; ...
    [-50, -45, alt_ref/1e3]; ...
    [ 50, -42, alt_ref/1e3]; ...
    [-60, -20, alt_ref/1e3];%[-90, -30, 140]; ...
    ]*1e3*scl.x;
v0 = [[1,0,0]; [-11,1,0]; [1,0,0]; [1,0,0]];%[0,0,1]
v1 = [[0,1,0]; [0,1,0]; [-2,3,0]; [0,1,0]];
r0 = [15, 20, 6, 1]*1e3*scl.x;
r1 = [3, 6, 11, 1]*1e3*scl.x;
colors = [ ...
    [1.0, 0.0, 0.0]; ...
    [0.0, 0.5, 0.0]; ...
    [1.0, 0.5, 0.0]; ...
    [1.0, 0.5, 0.7]; ...
    ];
res = [400,400,400,400]/2;
rev = [0,0,0,0];

views = [[-17-15,38];[90,0];[0,90]];
fluxes = zeros(2,ntubes);
v = 1;

% call flux tubes prior for speed
tubes = struct;
for n = 1:ntubes
    fprintf('Generating flux tube %i.\n',n)
    % tubes.(char(64+n)) = jules.tools.fluxtube(xg,dat,alt_ref*scl.x...
    %     ,p0(n,:),r0(n),r1(n),v0=v0(n,:),v1=v1(n,:)...
    %     ,reverse=rev(n),calculate_hull=0,res=res(n)...
    %     ,xlims = lim.x,ylims = lim.y);
    tubes.(char(64+n)) = jules.tools.current_flux_tube(xg,dat ...
        , p0(n,:), r0(n), r1(n), v0=v0(n,:), v1=v1(n,:) ...
        ,reverse=rev(n) ,res=res(n) ...
        ,xlims = lim.x, ylims = lim.y, zlims=lim.z);
end
vv = views(v,:);

fprintf('second: %i seconds.\n',round(toc))

figure(v)
set(gcf,'PaperUnits','inches','PaperPosition',[0,0,paper_w(v),paper_h(v)])
t = tiledlayout(1,1,'TileSpacing','compact');
axj = axes(t);
axn = axes(t);
axt = axes(t);
axa = [axj,axn,axt];

set(axa,'FontSize',fts)
set(axa,'FontName',ftn)
set(axj,'Color','none','XColor','none','YColor','none','ZColor','none','ZTick',alt_ref_p)
set(axn,'Color','none','XColor','none','YColor','none','ZColor','none','XTick',0)
set(axt,'Color','none','XGrid','on','YGrid','on','ZGrid','on')

xlim(axj,lim.x)
xlim(axn,lim.x + lim.x(1))
xlim(axt,lim.x)

ylim(axj,lim.y)
ylim(axn,lim.y)
ylim(axt,lim.y)

zlim(axj,lim.z + (alt_ref-zmin)*scl.x)
zlim(axn,lim.z)
zlim(axt,lim.z)

view(axa,vv)
ar = [range(lim.x),range(lim.y)*2,range(lim.z)];
pbaspect(axj,ar)
pbaspect(axn,ar)
pbaspect(axt,ar)
hold(axa,'on')

% fac slice plot
slice(axj,Xm_p,Ym_p,Zm_p,permute(jz_p,[2,1,3]),[],[],alt_ref_p);
colormap(axj,colorcet(clm.j))
shading(axj,'flat')
clim(axj,lim.j)
clb = colorbar(axj);
clb.Label.String = ['j_{||} (',units.j,')'];
clb.FontSize = fts;
clb.Position = [0.89,0.12,0.015,0.41];

% density slice plot
slice(axn,Xm_p,Ym_p,Zm_p,permute(log10(ne_p),[2,1,3]),0,[],[]);
colormap(axn,colorcet(clm.n))
clim(axn,[9.8,11.9])
shading(axn,'flat')
clb = colorbar(axn);
clb.Label.String = ['log n_e (',units.n,')'];
clb.FontSize = fts;
clb.Position = [0.89,0.57,0.015,0.41];
clb.Ruler.Exponent = clb_exp;

% plot flux tubes
for n = 1:ntubes
    % unpack tube
    color = colors(n,:);
    tube = tubes.(char(64+n));
    verts = tube.vertices;
    c0 = tube.caps.start;
    c1 = tube.caps.end;
    fluxes(1,n) = tube.flux.in*scl.j*scl.f;
    fluxes(2,n) = tube.flux.out*scl.j*scl.f;
    in0 = tube.flux.area.in;
    in1 = tube.flux.area.out;
    
    % generate top shadows
    shadow = nan(size(in0));
    shadow(in0) = 0;
    shadow(in1) = 0;
    shadow = repmat(shadow,[1,1,lz]);
    
    % generate east outlines
    pts = tube.points;
    outline_a = nan(2*lz-4,3);
    outline_b = nan(2*lz-4,3);
    for zid = 1:lz-2
        z0 = z(zid)*scl.x; z1 = z(zid+1)*scl.x;
        layer = pts(pts(:,3)>z0 & pts(:,3)<z1,:);        
        if size(layer,1) ~= 0
            [~,min_id] = min(layer(:,2));
            [~,max_id] = max(layer(:,2));
            outline_a(lz-zid-1,:) = layer(max_id,:);
            outline_a(lz-2+zid,:) = layer(min_id,:);

            layer_sorted = sort(layer',2)'; %#ok<UDIM>
            [gap,in_id] = max(diff(layer_sorted(:,2)));
            if gap > 1
                outline_b(lz-zid-1,:) = layer_sorted(in_id,:);
                outline_b(lz-2+zid,:) = layer_sorted(in_id+1,:);
            end
        end
    end
    outline_a = outline_a(not(isnan(outline_a(:,1))),:);
    outline_b = outline_b(not(isnan(outline_b(:,1))),:);

    pl0 = plot3(axt,c0(:,1),c0(:,2),c0(:,3));
    pl1 = plot3(axt,c1(:,1),c1(:,2),c1(:,3));
    stl = streamline(axt,verts);
    out_a = plot3(axn,zeros(size(outline_a(:,1))),outline_a(:,2),outline_a(:,3));
    out_b = plot3(axn,zeros(size(outline_b(:,1))),outline_b(:,2),outline_b(:,3));

    set(pl0,'Color','k','LineWidth',lnw*2)
    set(pl1,'Color','b','LineWidth',lnw*2)
    set(stl,'Color',[color,0.2],'LineWidth',lnw)
    set(out_a,'Color',color,'LineWidth',lnw)
    set(out_b,'Color',color,'LineWidth',lnw)
    if all(range(c0(:,3))<1)
        pl0s = plot3(axt,c0(:,1),c0(:,2),ones(size(c0(:,3)))*zmin*scl.x);
        set(pl0s,'Color','k','LineWidth',lnw*2,'LineStyle',':')
    end
    if all(range(c1(:,3))<1)
        pl1s = plot3(axt,c1(:,1),c1(:,2),ones(size(c1(:,3)))*zmin*scl.x);
        set(pl1s,'Color','b','LineWidth',lnw*2,'LineStyle',':')
    end
    if any(not(isnan(shadow)),'all')
        shd = slice(axj,Xm_p,Ym_p,Zm_p,permute(shadow,[2,1,3]),[],[],alt_ref_p);
        set(shd,'FaceAlpha',0.5)
    end
    shading(axj,'flat')
end

xlabel(sprintf('Mag. east (%s)',units.x))
ylabel(sprintf('Mag. north (%s)',units.x))
zlabel(sprintf('Mag. up (%s)',units.x))

flux_strings = cell(2,ntubes);
for n = 1:ntubes
    flux_strings{1,n} = ['\color[rgb]{',num2str(colors(n,:)),'}'...
        ,num2str(fluxes(1,n),'%3.1f'),' ',units.f,'  ','\color{black}'];
    flux_strings{2,n} = ['\color[rgb]{',num2str(colors(n,:)),'}'...
        ,num2str(fluxes(2,n),'%3.1f'),' ',units.f,'  ','\color{black}'];
end
if opts.annotate
    annotation('textbox',[0.52,0.97,0.01,0.01],'String',[...
        'Current In:  ',cell2mat(flux_strings(1,:)),newline,...
        'Current Out:  ',cell2mat(flux_strings(2,:)),newline,...
        ],'FitBoxToText','on','EdgeColor','none','BackgroundColor','none','FontSize',fts*18/20) %#ok<*UNRCH>
end

pan_y = 0.7; % >0 moves image down
zoom = opts.zoom;
for ax = axa
    campan(ax,pan_y,0);
    camzoom(ax,zoom);
end

fprintf('third: %i seconds.\n',round(toc))

filename = fullfile(pwd,'plots','paper0',['simulation-results',opts.suffix,'.png']);
saveas(gcf,filename)
fprintf('Saving: %s\n',filename)
% close all
end