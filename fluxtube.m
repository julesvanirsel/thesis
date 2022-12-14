% Generate a flux tube structure given a Gemini grid structure, Gemini data
% structure, and the central point, major, and minor radii of an ellipse.
% The output structure includes streamline vertices, a list of start and
% end cap coordinates, all coordinates defining the flux tube surface, and
% the area integrated value of the influx and outflux. Coordinates follow
% (x,y,z) = (East,North,Up).
%
% Example usage:
%   tube = fluxtube(xg,dat,300e3,[-1000e3,-50e3,300e3],200e3,20e3)
%
% Arguments:
%     xg                        Gemini grid structure
%     dat                       Gemini data structure
%     alt_ref                   reference altitude
%     p0                        ellipse centroid
%     r0                        ellipse major axis
%     r1                        ellipse minor axis
%     reverse                   (option) run flux tube in reverse
%     v0 = [1,0,0]              (option) ellipse major axis direction
%     v1 = [0,1,0]              (option) ellipse minor axis direction
%     plot = false              (option) plotting option
%     show_reverse_test = false (option) plot reverse flux tube test
%     calculate_hull = false    (option) calculate flux tube hull
%     var = 'j'                 (option) flux tube variable ('j','v','b')
%     res = 32                  (option) number of flux tube streamlines
%     color = [0,0,0]           (option) plotting color
%     units = 'km'              (option) length units ('m','km')
%
% Dependencies:
%   matlab R2020a or higher
%
% Contact: jules.van.irsel.gr@dartmouth.edu

function tube = fluxtube(xg,dat,alt_ref,p0,r0,r1,options)
arguments
    xg (1,1) struct {mustBeNonempty}
    dat (1,1) struct {mustBeNonempty}
    alt_ref (1,1) double {mustBePositive}
    p0 (1,3) double {mustBeNonempty}
    r0 (1,1) double {mustBeNonempty}
    r1 (1,1) double {mustBeNonempty}
    options.reverse (1,1) logical {mustBeNonempty} = false
    options.v0 (1,3) double {mustBeNonempty} = [1,0,0]
    options.v1 (1,3) double {mustBeNonempty} = [0,1,0]
    options.xlims (1,2) double {mustBeNonempty} = [-1,1]*1500;
    options.ylims (1,2) double {mustBeNonempty} = [-1,1]*400;
    options.zmin (1,1) double {mustBeNonnegative} = 80;
    options.plot (1,1) logical {mustBeNonempty} = false
    options.show_reverse_test (1,1) logical {mustBeNonempty} = false
    options.show_flux (1,1) logical {mustBeNonempty} = false
    options.calculate_hull (1,1) logical {mustBeNonempty} = false
    options.var (1,1) char {mustBeMember(options.var,['j','v','b'])} = 'j'
    options.res (1,1) int16 {mustBePositive} = 32
    options.color (1,3) double {mustBeNonempty} = [0,0,0]
    options.units (1,:) char {mustBeMember(options.units,['m','km'])} = 'km'
end
%% preamble
if strcmp(options.units,'km')
    x_scl = 1e-3;
elseif strcmp(options.units,'m')
    x_scl = 1;
end

%% define coordinates
x = double(xg.x2(3:end-2)*x_scl);
y = double(xg.x3(3:end-2)*x_scl);
z = double(xg.x1(3:end-2)*x_scl);
[~,ubz] = min(abs(z-alt_ref));
ubz = ubz + 1; % add buffer cell
z = z(1:ubz);
dx = double(xg.dx2h*x_scl);
dy = double(xg.dx3h*x_scl);
dz = double(xg.dx1h(1:ubz)*x_scl);
dmin = min([dx;dy;dz]);
lx = length(x); ly = length(y); lz = length(z);
xlims = options.xlims;
ylims = options.ylims;
zmin = options.zmin;
[I,J,K] = ndgrid(1:lx,1:ly,1:lz);
[X,Y,Z] = ndgrid(x,y,z);
[Xm,Ym,Zm] = meshgrid(x,y,z);
[dX,dY,~] = ndgrid(dx,dy,dz);
dA = dX.*dY;
gi = griddedInterpolant(X,Y,Z,I); % stream3 is more robust in index space. back and forth interpolation required
gj = griddedInterpolant(X,Y,Z,J);
gk = griddedInterpolant(X,Y,Z,K);
gx = griddedInterpolant(I,J,K,X);
gy = griddedInterpolant(I,J,K,Y);
gz = griddedInterpolant(I,J,K,Z);

%% calculate initial curve in index space
v0 = options.v0/norm(options.v0);
v1 = options.v1/norm(options.v1);
t = repmat(linspace(0,1,options.res)',1,3);
c0 = p0 + r0*cos(2*pi*t).*v0 + r1*sin(2*pi*t).*v1;
c0_ids = [gi(c0),gj(c0),gk(c0)];

%% assert all initial coordinates are inside the model space
assert(all(1 <= c0_ids(:,1)) && all(c0_ids(:,1) <= lx),'Initial flux tube curve out of bounds.')
assert(all(1 <= c0_ids(:,2)) && all(c0_ids(:,2) <= ly),'Initial flux tube curve out of bounds.')
assert(all(1 <= c0_ids(:,3)) && all(c0_ids(:,3) <= lz),'Initial flux tube curve out of bounds.')

%% gather flux tube vector field
if strcmp(options.var,'j')
    Vx = permute(dat.J2,[3,2,1]); % stream3 requires meshgrid format
    Vy = permute(dat.J3,[3,2,1]);
    Vz = permute(dat.J1,[3,2,1]);
elseif strcmp(options.var,'v')
    Vx = permute(dat.v2,[3,2,1]); % stream3 requires meshgrid format
    Vy = permute(dat.v3,[3,2,1]);
    Vz = permute(dat.v1,[3,2,1]);
elseif strcmp(options.var,'b')
    error('Magnetic field fluxtubes not yet incorperated')
end

%% change vector field boundaries
Vx = Vx(:,:,1:ubz);
Vy = Vy(:,:,1:ubz);
Vz = Vz(:,:,1:ubz);

%% calculate flux tube streamline vertices
if options.reverse
    dir = -1;
else
    dir = 1;
end
step = 1e-2;
maxvert = 1e6;
verts_ids = stream3(dir*Vx,dir*Vy,dir*Vz,c0_ids(:,1),c0_ids(:,2),c0_ids(:,3),[step maxvert]);
verts = cellfun(@(v) [gx(v),gy(v),gz(v)],verts_ids,'UniformOutput',0);
pts = cell2mat(cellfun(@(v) v',verts,'UniformOutput',0))';
c1 = cell2mat(cellfun(@(v) v(end,:)',verts,'UniformOutput',0))';

% matching reverse flux tube test
c1_ids = [gi(c1),gj(c1),gk(c1)];
verts_ids_test = stream3(-dir*Vx,-dir*Vy,-dir*Vz,c1_ids(:,1),c1_ids(:,2),c1_ids(:,3),[step maxvert]);
verts_test = cellfun(@(v) [gx(v),gy(v),gz(v)],verts_ids_test,'UniformOutput',0);
c0_test = cell2mat(cellfun(@(v) v(end,:)',verts_test,'UniformOutput',0))';
if any(abs(c0(:,1:2)-c0_test(:,1:2)) > dmin)
    warning('Reverse flux tube test failed.')
end

%% calculate in- and outflux (limited to horizontal endcaps)
c0_cid = round(mean(c0_ids));
c1_cid = round(mean(c1_ids));
in0 = inpolygon(squeeze(X(:,:,c0_cid(3))),squeeze(Y(:,:,c0_cid(3))),c0(:,1),c0(:,2));
in1 = inpolygon(squeeze(X(:,:,c1_cid(3))),squeeze(Y(:,:,c1_cid(3))),c1(:,1),c1(:,2));
flux0 = sum(squeeze(Vz(:,:,c0_cid(3)))'.*squeeze(-dir*dA(:,:,c0_cid(3))).*in0,'all');
flux1 = sum(squeeze(Vz(:,:,c1_cid(3)))'.*squeeze(dir*dA(:,:,c1_cid(3))).*in1,'all');
if max(c0(:,3))-min(c0(:,3)) > dmin
    warning('Nonhorizontal start curve may result in inaccurate influx value.')
end
if max(c1(:,3))-min(c1(:,3)) > dmin
    warning('Nonhorizontal end curve may result in inaccurate outflux value.')
end

%% calculate flux tube hull
if options.calculate_hull
    hull = false(size(I));
    for h = 1:lz
        ptsA = [];
        ptsB = [];
        for n = 1:length(verts_ids)
            v = verts_ids{n};
            [~,mid] = min(v(:,3)); % split streamline by min. height vertex
            vA = v(1:mid,:);
            vB = v(mid+1:end,:);
            ptsA = [ptsA;vA(round(vA(:,3))==h,1:2)]; %#ok<AGROW> % collect points around height h
            ptsB = [ptsB;vB(round(vB(:,3))==h,1:2)]; %#ok<AGROW>
        end
        if ~isempty(ptsA) && ~isempty(ptsB)
            ptsA = ptsA(convhull(ptsA),:); % grab 2d convex hull of points at height h
            ptsB = ptsB(convhull(ptsB),:);
        end
        inA = inpolygon(squeeze(I(:,:,h)),squeeze(J(:,:,h)),ptsA(:,1),ptsA(:,2));
        inB = inpolygon(squeeze(I(:,:,h)),squeeze(J(:,:,h)),ptsB(:,1),ptsB(:,2));
        hull(:,:,h) = inA | inB;
    end
end

%% plotting
if options.plot
    sized = 20;
    zlims = [zmin,alt_ref*1.05]; % add 5% buffer
    shadow = nan(size(in0));
    shadow(in0) = 0;
    shadow(in1) = 0;

    figure(1)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    hold on
    pl0 = plot3(c0(:,1),c0(:,2),c0(:,3));
    pl1 = plot3(c1(:,1),c1(:,2),c1(:,3));
    stl = streamline(verts);
    pcolor(squeeze(X(:,:,end)),squeeze(Y(:,:,end)),squeeze(Vz(:,:,end)'))
    shd = pcolor(squeeze(X(:,:,end)),squeeze(Y(:,:,end)),shadow);
    if options.calculate_hull
        isosurface(Xm,Ym,Zm,permute(hull*mean(Vz(:)),[2,1,3]),mean(Vz(:)))
    end
    if options.show_reverse_test
        pl0_test = plot3(c0_test(:,1),c0_test(:,2),c0_test(:,3));
        stl_test = streamline(verts_test);
        set(pl0_test,'Color','r','LineWidth',2)
        set(stl_test,'Color','r','LineWidth',0.5)
    end
    hold off
    set(shd,'FaceAlpha',0.5)
    shading flat
    clb = colorbar;
    clb.Label.String = ['Vz (',num2str(alt_ref),' ',options.units,')'];
    clb.FontSize = sized;
    set(pl0,'Color',options.color,'LineWidth',2)
    set(pl1,'Color',options.color,'LineWidth',2)
    set(stl,'Color',options.color,'LineWidth',0.5)
    set(gca,'FontSize',sized)
    view([30,45])
    xlim(xlims)
    ylim(ylims)
    zlim(zlims)
    ar = [0.4*range(xlims),range(ylims),2*range(zlims)]; %%%CHANGE
    %     ar = [0.2*range(xlims),range(ylims),range(zlims)]; %%%CHANGE
    pbaspect(ar)
    xlabel(['East [',options.units,']'],'FontSize',sized)
    ylabel(['North [',options.units,']'],'FontSize',sized)
    zlabel(['Up [',options.units,']'],'FontSize',sized)
    if options.show_flux
        text(c0(1,1)-6e2,c0(1,2)+0e2,1.05*alt_ref,[num2str(flux0,3),'<V> ',options.units,'^2'],'Color',options.color,'FontSize',0.8*sized)
        text(c1(1,1)+1e2,c1(1,2)+1e2,1.05*alt_ref,[num2str(flux1,3),'<V> ',options.units,'^2'],'Color',options.color,'FontSize',0.8*sized)
    end
end

tube.vertices = verts;
tube.verts_ids = verts_ids;
tube.caps.start = c0;
tube.caps.end = c1;
if options.calculate_hull
    tube.hull = hull;
end
tube.points = pts;
tube.flux.in = flux0;
tube.flux.out = flux1;
tube.flux.area.in = in0;
tube.flux.area.out = in1;
end