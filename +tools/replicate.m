% Description:
%   aaa
%
% Example usage:
%   phi = replicator()
%
% Arguments:
%   xdata
%
% Dependencies:
%   matlab R2022a or higher
%   gemini3d
%   colorcet
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu

function [phi,v2_int,v3_int] = replicate(in_situ,image,xg,opts)
arguments
    in_situ (1,1) struct {mustBeNonempty}
    image (1,1) struct {mustBeNonempty}
    xg (1,1) struct {mustBeNonempty}
    opts.num_replications (1,1) int32 {mustBePositive} = 256
    opts.pos_type {mustBeMember(opts.pos_type,["angular","linear"])} = "angular"
    opts.flow_bg (1,2) double {mustBeNumeric} = [nan,nan]
    opts.flow_smoothing_window (1,1) int32 {mustBePositive} = 8
    opts.show_plots (1,1) logical = false
    opts.save_plots (1,1) logical = false
    opts.do_rotate (1,1) logical = true
    opts.do_scale (1,1) logical = true
    opts.auto_lim (1,1) logical = true
    opts.direc (1,:) char {mustBeFolder} = '.'
end

%% initialize
ftn = 'Arial';
fts = 11;

close all
reset(0)
set(0,'defaultFigurePaperUnits','inches')
set(0,'defaultTiledlayoutTileSpacing','tight')
set(0,'defaultSurfaceEdgeColor','flat')
tools.setall(0,'FontName',ftn)
tools.setall(0,'FontSize',fts)
tools.setall(0,'Multiplier',1)

scl.x = 1e-3;   unt.x = 'km';       clm.dv = 'CBD1';
scl.v = 1e-3;   unt.v = 'km/s';     clm.v = 'D2';
scl.Q = 1e0;    unt.Q = 'mW/m^2';   clm.Q = 'L19';
scl.p = 1e-3;   unt.p = 'kV';       clm.p = 'D10';

lim.y = [-1,1]*59;  lim.x = [-1,1]*120;
lim.v = [-1,1]*1.3; lim.dv = [-1,1]*0.1;

lbl.x = sprintf('Magnetic east (%s)',unt.x);
lbl.y = sprintf('Magnetic north (%s)',unt.x);
lbl.vx = sprintf('Eastward flow (%s)',unt.v);
lbl.vy = sprintf('Northward flow (%s)',unt.v);

% assertions
assert(isequal(size(in_situ.pos),size(in_situ.flow)), ...
    'in_situ.pos and in_situ.flow must have equal sizes.')
assert(isequal(size(image.pos,1:2),size(image.flux)), ...
    'image.flux and image.pos(:,:,1) must have equal sizes.')

% unpack grid
MLAT = 90-squeeze(xg.theta(end,:,:))*180/pi;
MLON = squeeze(xg.phi(end,:,:))*180/pi;
x2 = double(xg.x2(3:end-2))';
x3 = double(xg.x3(3:end-2))';
dx2 = xg.dx2h;
dx3 = xg.dx3h;
lx2 = xg.lx(2); lx3 = xg.lx(3);
[X2,X3] = ndgrid(x2,x3);
[DX2,DX3] = ndgrid(dx2,dx3);
mlon_to_x2 = griddedInterpolant(MLON(:,1),x2);
mlat_to_x3 = griddedInterpolant(MLAT(1,:),x3);
Bmag = abs(mean(xg.Bmag,'all'));

% unpack in situ data + image
if strcmp(opts.pos_type,"angular")
    mlon_traj = in_situ.pos(:,1);
    mlat_traj = in_situ.pos(:,2);
    x2_traj = mlon_to_x2(mlon_traj);
    x3_traj = mlat_to_x3(mlat_traj);

    mlon_imag = image.pos(:,:,1);
    mlat_imag = image.pos(:,:,2);
    X2_imag = mlon_to_x2(mlon_imag);
    X3_imag = mlat_to_x3(mlat_imag);
else
    x2_traj = in_situ.pos(:,1);
    x3_traj = in_situ.pos(:,2);
    X2_imag = image.pos(:,:,1);
    X3_imag = image.pos(:,:,2);
end

v2_traj = in_situ.flow(:,1);
v3_traj = in_situ.flow(:,2);
[~,sort_ids] = sort(x2_traj);

Q = image.flux;
x2_imag = X2_imag(:,1)';
x3_imag = X3_imag(1,:);

%% calculate boundaries
edges = tools.find_max_edges(Q);
num_bounds = 2;
bounds = nan(num_bounds,size(edges,1));
for i = 1:size(edges,1)
    [~,bounds(:,i)] = tools.peak_detect(edges(i,:), ...
        num=num_bounds,smoothness=0.009);
end
x2_bounds = sort(x2_imag(2:end-1));
x3_bounds = smoothdata(x3_imag(bounds+1),2,'gaussian');
bound.A = griddedInterpolant(x2_bounds,x3_bounds(1,:));
bound.B = griddedInterpolant(x2_bounds,x3_bounds(2,:));
angle = griddedInterpolant(x2_bounds(1:end-1), ...
    atan2(diff(x3_bounds(1,:)),diff(x2_bounds)));
bound_pts = linspace(min(x2_imag),max(x2_imag),length(x2_imag));

if opts.show_plots
    for p = 1:3
        figure(p)
        hold on
        if p==1
            pcolor(X2_imag*scl.x,X3_imag*scl.x,Q*scl.Q)
            clb = colorbar;
            clb.Label.String = sprintf('Total energy flux(%s)',unt.Q);
        elseif p==2
            contour(X2_imag*scl.x,X3_imag*scl.x,Q*scl.Q)
        else
            pcolor(X2_imag(2:end-1,2:end-1)*scl.x,X3_imag(2:end-1,2:end-1)*scl.x,edges)
        end
        scatter(x2_bounds*scl.x,x3_bounds(1,:)*scl.x,'k')
        scatter(x2_bounds*scl.x,x3_bounds(2,:)*scl.x,'b')
        plot(bound_pts*scl.x,bound.A(bound_pts)*scl.x,'k')
        plot(bound_pts*scl.x,bound.B(bound_pts)*scl.x,'b')
        colormap(colorcet(clm.Q))
        xlabel(lbl.x); ylabel(lbl.y)
        xlim(lim.x); ylim(lim.y)

        legend('','Boundary A','Boundary B')
    end
end

%% replicate in situ flow data
if all(isnan(opts.flow_bg))
    v_bg = mean(in_situ.flow);
else
    v_bg = opts.flow_bg;
end
v2_traj = smoothdata(v2_traj - v_bg(1),'gaussian',opts.flow_smoothing_window);
v3_traj = smoothdata(v3_traj - v_bg(2),'gaussian',opts.flow_smoothing_window);

% 0 = original, 1 = replicated, a = primary bound, b = secondary bound
traj0 = griddedInterpolant(tools.minsmooth(x2_traj(sort_ids)),x3_traj(sort_ids)); % smooth needed for unique
x0a = fzero(@(x)(traj0(x)-bound.A(x)),0);
x0b = fzero(@(x)(traj0(x)-bound.B(x)),0);
y0a = traj0(x0a);
y0b = traj0(x0b);
% y0a = bound.A(mean([x0a,x0b]));
% y0b = bound.B(mean([x0a,x0b]));
width0 = sqrt((x0b-x0a)^2+(y0b-y0a)^2);
% width0 = abs(y0b-y0a);
beta = atan2(x0b-x0a,y0b-y0a); % angle b/w bound-traj inters. and vertical

dx_min = (x2(1) - max(x2_traj))*1.1;
dx_max = (x2(end) - min(x2_traj))*1.1;
dxs = linspace(dx_min,dx_max,opts.num_replications); % eastward displacements
x2_traj_rep = nan(length(dxs),length(x2_traj));
x3_traj_rep = nan(length(dxs),length(x3_traj));
v2_traj_rep = nan(length(dxs),length(v2_traj));
v3_traj_rep = nan(length(dxs),length(v3_traj));
for i = 1:length(dxs)
    % translate
    dx = dxs(i);
    dy = bound.A(x0a+dx)-y0a;
    x2_traj_tra = x2_traj + dx;
    x3_traj_tra = x3_traj + dy;

    % determine width at position 1
    traj1 = griddedInterpolant(tools.minsmooth(x2_traj_tra(sort_ids)),x3_traj_tra(sort_ids));
    x1a = fzero(@(x)(traj1(x)-bound.A(x)),0);
    x1b = fzero(@(x)(traj1(x)-bound.B(x)),0);
    y1a = traj1(x1a);
    y1b = traj1(x1b);
%     y1a = bound.A(mean([x1a,x1b]));
%     y1b = bound.B(mean([x1a,x1b]));
    width1 = sqrt((x1b-x1a)^2+(y1b-y1a)^2);
%     width1 = abs(y1b-y1a);

    if opts.do_scale
%         beta = angle(x1a);

        % rotate about p1a by beta
        x2_traj_rot = cos(beta)*(x2_traj_tra - x1a) ...
            - sin(beta)*(x3_traj_tra - y1a) + x1a;
        x3_traj_rot = sin(beta)*(x2_traj_tra - x1a) ...
            + cos(beta)*(x3_traj_tra - y1a) + y1a;

        % scale about p1a
        scale = width1/width0;
        x2_traj_scl = x2_traj_rot;
        x3_traj_scl = scale*(x3_traj_rot - y1a) + y1a;

        % rotate back
        x2_traj_rep(i,:) = cos(-beta)*(x2_traj_scl - x1a) ...
            - sin(-beta)*(x3_traj_scl - y1a) + x1a;
        x3_traj_rep(i,:) = sin(-beta)*(x2_traj_scl - x1a) ...
            + cos(-beta)*(x3_traj_scl - y1a) + y1a;
    else
        x2_traj_rep(i,:) = x2_traj_tra;
        x3_traj_rep(i,:) = x3_traj_tra;
    end

    if opts.do_rotate
        % rotate flows to be tangent to bound_a
        chi = angle(x1a);
        v2_traj_rep(i,:) = cos(chi)*v2_traj - sin(chi)*v3_traj;
        v3_traj_rep(i,:) = sin(chi)*v2_traj + cos(chi)*v3_traj;
    else
        v2_traj_rep(i,:) = v2_traj;
        v3_traj_rep(i,:) = v3_traj;
    end

    if i==50 && opts.show_plots
        figure(99)
        hold on
        pcolor(X2_imag*scl.x,X3_imag*scl.x,Q*scl.Q)
        plot(bound_pts*scl.x,bound.A(bound_pts)*scl.x,'b')
        plot(bound_pts*scl.x,bound.B(bound_pts)*scl.x,'b')
        quiver(x2_traj_rot'*scl.x,x3_traj_rot'*scl.x, ...
            v2_traj_rep(i,:)*scl.v,v3_traj_rep(i,:)*scl.v,1,'.-k')
        quiver(x2_traj_rep(i,:)*scl.x,x3_traj_rep(i,:)*scl.x, ...
            v2_traj_rep(i,:)*scl.v,v3_traj_rep(i,:)*scl.v,1,'.-k')
        quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,1,'.-g')
        scatter([x0a,x0b,x1a,x1b]*scl.x,[y0a,y0b,y1a,y1b]*scl.x,'Filled','k')
        xlim(lim.x); ylim(lim.y)
        xlabel(lbl.x); ylabel(lbl.y)
        colormap(colorcet(clm.Q))
        clb = colorbar;
        clb.Label.String = sprintf('Total energy flux(%s)',unt.Q);
    end
end

if opts.show_plots
    figure(4)
    hold on
    plot(x3_traj*scl.x,v2_traj*scl.v,'r')
    plot(x3_traj*scl.x,v3_traj*scl.v,'b')
    plot(x3_traj*scl.x,(in_situ.flow(:,1) - v_bg(1))*scl.v,'Color',[1,0.4,0.4])
    plot(x3_traj*scl.x,(in_situ.flow(:,2) - v_bg(2))*scl.v,'Color',[0.4,0.4,1])
    xlabel(lbl.y); ylabel(sprintf('Flow (%s)',unt.v))
    legend('Eastward','Eastward smoothed','Northward','Northward smoothed')

    figure(5)
    hold on
    pcolor(X2_imag*scl.x,X3_imag*scl.x,Q*scl.Q)
    plot(bound_pts*scl.x,bound.A(bound_pts)*scl.x,'b')
    plot(bound_pts*scl.x,bound.B(bound_pts)*scl.x,'b')
    quiver(x2_traj_rep*scl.x,x3_traj_rep*scl.x, ...
        v2_traj_rep*scl.v,v3_traj_rep*scl.v,0,'.-k')
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,0,'.-g')
    xlim(lim.x); ylim(lim.y)
    xlabel(lbl.x); ylabel(lbl.y)
    colormap(colorcet(clm.Q))
    clb = colorbar;
    clb.Label.String = sprintf('Total energy flux(%s)',unt.Q);
end

%% interpolate replicated flow data
fv2 = scatteredInterpolant(x2_traj_rep(:),x3_traj_rep(:),v2_traj_rep(:));
fv3 = scatteredInterpolant(x2_traj_rep(:),x3_traj_rep(:),v3_traj_rep(:));
v2_int = fv2(X2,X3);
v3_int = fv3(X2,X3);
E2_int =  v3_int*Bmag; % v = ExB/B^2
E3_int = -v2_int*Bmag;

%% generate potential map
% translated from Alex Mule's python code Oct 9, 2023

% extrapolate data to avoid edge effects with fourier transforms
nfill = 16;
x2_ext = [x2(1)+dx2(1)*(-nfill:-1),x2,x2(end)+dx2(end)*(1:nfill)];
x3_ext = [x3(1)+dx3(1)*(-nfill:-1),x3,x3(end)+dx3(end)*(1:nfill)];
[X2_ext,X3_ext] = ndgrid(x2_ext,x3_ext);
fE2_int = griddedInterpolant(X2,X3,E2_int,'linear','nearest');
fE3_int = griddedInterpolant(X2,X3,E3_int,'linear','nearest');
E2_int_ext = fE2_int(X2_ext,X3_ext);
E3_int_ext = fE3_int(X2_ext,X3_ext);

% wave vector convention: 2 pi ( -f_Ny : +f_Ny )
k2_Ny = 2*pi/(2*mean(dx2));
k3_Ny = 2*pi/(2*mean(dx3));
k2 = linspace(-k2_Ny,k2_Ny,lx2+2*nfill);
k3 = linspace(-k3_Ny,k3_Ny,lx3+2*nfill);
[K2,K3] = ndgrid(k2,k3);

% Fourier transform of electric field
G2 = fftshift(fft2(E2_int_ext));
G3 = fftshift(fft2(E3_int_ext));

% Fourier transform of phi
Gphi = 1i*(K2.*G2+K3.*G3)./(K2.^2+K3.^2);

% inverse Fourier transform of Gphi + resampling onto working grid
phi0 = real(ifft2(ifftshift(Gphi)));
phi0 = phi0(nfill+1:end-nfill,nfill+1:end-nfill);

% determine average electric field of phi0
if range(dx2) < 1e-3
    [E20,E30] = gradient(-phi0',mean(dx2),mean(dx3));
else
    [E20,E30] = gradient(-phi0',dx2,dx3);
end
E20 = E20';
E30 = E30';

% make average electric field match
phi = phi0 - mean(E2_int-E20,'all').*X2 - mean(E3_int-E30,'all').*X3;
phi = phi - mean(phi,'all');

if opts.show_plots
    figure(6)
    pcolor(X2*scl.x,X3*scl.x,phi*scl.p)
    colormap(colorcet(clm.p))
    clb = colorbar;
    clb.Label.String = sprintf('Electric potential (%s)',unt.p);
end

%% plot final flow fields
if opts.show_plots || opts.save_plots
    if range(dx2) < 1e-3
        [E2,E3] = gradient(-phi',mean(dx2),mean(dx3));
    else
        [E2,E3] = gradient(-phi',dx2,dx3);
    end
    v2 = -E3'/Bmag;
    v3 =  E2'/Bmag;
    divv_int = diff(v2_int(:,2:end),1,1)./DX2(2:end,2:end) + ...
        diff(v3_int(2:end,:),1,2)./DX3(2:end,2:end);
    divv = diff(v2(:,2:end),1,1)./DX2(2:end,2:end) + ...
        diff(v3(2:end,:),1,2)./DX3(2:end,2:end);

    figure(7)
    set(gcf,'PaperPosition',[0,0,9,6.5])
    tiledlayout(3,3);

    if opts.auto_lim
        qnt = 0.95;
        max_v = quantile(abs([v2(:)+v_bg(1);v3(:)+v_bg(2)]),qnt);
        max_dv = quantile(abs([divv(:);divv_int(:)]),qnt);
        max_p = quantile(abs(phi(:)),qnt);
        lim.v = [-1,1]*max_v*scl.v;
        lim.dv = [-1,1]*max_dv;
        lim.p = [-1,1]*max_p;
    end

    % row 1
    nexttile
    title('Interpolated')
    hold on
    pcolor(X2*scl.x,X3*scl.x,(v2_int+v_bg(1))*scl.v)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
    colormap(gca,colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = sprintf('Eastward flow (%s)',unt.v);
    clb.Location = 'westoutside';
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([])
    ylabel(lbl.y)

    nexttile
    title('Helmholtz decomposition')
    hold on
    pcolor(X2*scl.x,X3*scl.x,(v2+v_bg(1))*scl.v)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
    colormap(gca,colorcet(clm.v))
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([]); yticks([])

    nexttile
    title('Difference & Potential map')
    hold on
    pcolor(X2*scl.x,X3*scl.x,(v2-v2_int)*scl.v)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
    colormap(gca,colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = sprintf('\\Delta Eastward flow (%s)',unt.v);
    xlim(lim.x); ylim(lim.y); clim(lim.v/3)
    xticks([]); yticks([])

    %row 2
    nexttile
    hold on
    pcolor(X2*scl.x,X3*scl.x,(v3_int+v_bg(2))*scl.v)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
    colormap(gca,colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = sprintf('Northward flow (%s)',unt.v);
    clb.Location = 'westoutside';
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([])
    ylabel(lbl.y)

    nexttile
    hold on
    pcolor(X2*scl.x,X3*scl.x,(v3+v_bg(2))*scl.v)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
    colormap(gca,colorcet(clm.v))
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([]); yticks([])

    nexttile
    hold on
    pcolor(X2*scl.x,X3*scl.x,(v3-v3_int)*scl.v)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
    colormap(gca,colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = sprintf('\\Delta Northward flow (%s)',unt.v);
    xlim(lim.x); ylim(lim.y); clim(lim.v/3)
    xticks([]); yticks([])

    % row 3
    nexttile
    hold on
    pcolor(X2(1:end-1,1:end-1)*scl.x,X3(1:end-1,1:end-1)*scl.x,divv_int)
    colormap(gca,colorcet(clm.dv))
    clb = colorbar;
    clb.Label.String = 'Divergence (Hz)';
    clb.Location = 'westoutside';
    xlim(lim.x); ylim(lim.y); clim(lim.dv)
    xlabel(lbl.x); ylabel(lbl.y)

    nexttile
    hold on
    pcolor(X2(1:end-1,1:end-1)*scl.x,X3(1:end-1,1:end-1)*scl.x,divv)
    colormap(gca,colorcet(clm.dv))
    xlim(lim.x); ylim(lim.y); clim(lim.dv)
    yticks([])
    xlabel(lbl.x)

    nexttile
    hold on
    pcolor(X2*scl.x,X3*scl.x,phi*scl.p)
    colormap(gca,colorcet(clm.p))
    clb = colorbar;
    clb.Label.String = sprintf('Electric potential (%s)',unt.p);
    xlim(lim.x); ylim(lim.y);
    yticks([])
    xlabel(lbl.x)
    
    if opts.save_plots
        filename = sprintf('replicated_%s.png',datetime(datetime,'Format','MMM-dd''_''hh-mm-ss'));
        saveas(gcf,fullfile(opts.direc,filename))
    end
    if ~opts.show_plots
        close all
    end
end
end