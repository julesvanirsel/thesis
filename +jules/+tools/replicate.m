% Description:
%   aaa
%
% Example usage:
%   phi = replicate()
%
% Arguments:
%   in_situ
%   image
%   xg
%
% Dependencies:
%   matlab R2022a or higher
%   gemini3d
%   colorcet
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu

function [phi,mlon,mlat,E2_bg,E3_bg,v2_int,v3_int] = replicate(in_situ,image,xg,opts)
arguments
    in_situ (1,1) struct {mustBeNonempty}
    image (1,1) struct {mustBeNonempty}
    xg (1,1) struct {mustBeNonempty}
    opts.upsample (1,1) int16 {mustBePositive} = 1
    opts.num_replications (1,1) int32 {mustBePositive} = 256
    opts.pos_type {mustBeMember(opts.pos_type,["angular","linear"])} = "angular"
    opts.flow_bg (1,2) double {mustBeNumeric} = [nan,nan]
    opts.flow_smoothing_window (1,1) int32 {mustBePositive} = 1
    opts.do_rotate (1,1) logical = true
    opts.do_scale (1,1) logical = true
    opts.arc_definition {mustBeMember(opts.arc_definition,["conductance","flux"])} = "conductance"
    opts.edge_method {mustBeMember(opts.edge_method,["contour","sobel"])} = "contour"
    opts.contour_values (1,2) double = [nan,nan]
    opts.boundary_smoothing_window (1,1) int32 {mustBePositive} = 1
    opts.swap_primary (1,1) logical = false
    opts.fit_harmonic (1,1) logical = true
    opts.harmonic_mask (1,3) double {mustBeNonnegative} = [1,1,2]*5e3
    opts.add_phi_background (1,1) logical = false
    opts.plot_bg (1,1) logical = false
    opts.show_plots (1,1) logical = false
    opts.save_plots (1,3) logical = false(3,1)
    opts.save_data (1,1) logical = false
    opts.auto_lim (1,1) logical = true
    opts.direc (1,:) char {mustBeFolder} = 'plots'
    opts.suffix (1,:) char = ''
    opts.starting_letter (1,1) char = 'A'
end

%% initialize
ftn = 'Arial';
fts = 10*2;
lw = 1.4;

close all
reset(0)
set(0,'defaultFigurePaperUnits','inches')
set(0,'defaultTiledlayoutPadding','tight')
set(0,'defaultTiledlayoutTileSpacing','tight')
set(0,'defaultSurfaceEdgeColor','flat')
set(0,'defaultLineLineWidth',lw)
set(0,'defaultScatterLineWidth',lw)
set(0,'defaultQuiverLineWidth',lw*0.7)
jules.tools.setall(0,'FontName',ftn)
jules.tools.setall(0,'FontSize',fts)
jules.tools.setall(0,'Multiplier',1)

colorcet = @jules.tools.colorcet;

scl.x = 1e-3; scl.v = 1e-3; scl.dv = 1e3; scl.p = 1e-3;
unt.x = 'km'; unt.v = 'km/s'; unt.dv = 'mHz'; unt.p = 'kV';
clm.v = 'D2'; clm.dv = 'CBD1'; clm.p = 'D10';
lim.x = [-1,1]*125; lim.y = [-1,1]*59;  lim.v = [-1,1]*1.3; lim.dv = [-1,1]*0.1;

lbl.x = sprintf('Mag. E (%s)',unt.x);
lbl.y = sprintf('Mag. N (%s)',unt.x);
lbl.vx = sprintf('v_E (%s)',unt.v);
lbl.vy = sprintf('v_N (%s)',unt.v);

if not(isempty(opts.suffix))
    opts.suffix = ['_',opts.suffix];
end

% assertions
assert(isequal(size(in_situ.pos),size(in_situ.flow)), ...
    'in_situ.pos and in_situ.flow must have equal sizes.')
assert(isequal(size(image.pos,1:2),size(image.flux)), ...
    'image.flux and image.pos(:,:,1) must have equal sizes.')
assert(size(in_situ.pos,2)==2, ...
    'Second dimension of in_situ.pos and in_situ.flow must be 2.')
assert(size(image.pos,3)==2, ...
    'Third dimension of image.pos must be 2.')
assert(issorted(opts.contour_values), ...
    'Contour values must be in ascending order.')

% unpack grid
if opts.upsample > 1
    n_us = opts.upsample;
    fprintf('Upsampling grid by a factor of %i.\n',n_us)
    theta_tmp = linspace(max(xg.theta(1,1,:)),min(xg.theta(1,1,:)),xg.lx(3)*n_us);
    phi_tmp = linspace(min(xg.phi(1,:,1)),max(xg.phi(1,:,1)),xg.lx(2)*n_us);
    xg_us.theta = permute(repmat(theta_tmp,[xg.lx(1),1,xg.lx(2)*n_us]),[1,3,2]);
    xg_us.phi = permute(repmat(phi_tmp,[xg.lx(1),1,xg.lx(3)*n_us]),[1,2,3]);
    xg_us.x2 = linspace(min(xg.x2),max(xg.x2),(length(xg.x2)-4)*n_us+4);
    xg_us.x3 = linspace(min(xg.x3),max(xg.x3),(length(xg.x3)-4)*n_us+4);
    xg_us.dx2h = ones(1,numel(xg.dx2h)*n_us)*mean(xg.dx2h);
    xg_us.dx3h = ones(1,numel(xg.dx3h)*n_us)*mean(xg.dx3h);
    xg_us.lx = [xg.lx(1), xg.lx(2)*n_us, xg.lx(3)*n_us];
    xg_us.Bmag = xg.Bmag;
    xg = xg_us;
end

mlat = 90-squeeze(xg.theta(1,1,:))*180/pi;
mlon = squeeze(xg.phi(1,:,1))*180/pi;
x2 = double(xg.x2(3:end-2));
x3 = double(xg.x3(3:end-2));
if size(x2,1)~=1; x2=x2'; end
if size(x3,1)~=1; x3=x3'; end
dx2 = xg.dx2h;
dx3 = xg.dx3h;
lx2 = xg.lx(2); lx3 = xg.lx(3);
[X2,X3] = ndgrid(x2,x3);
mlon_to_x2 = griddedInterpolant(mlon,x2);
mlat_to_x3 = griddedInterpolant(mlat,x3);
Bmag = abs(mean(xg.Bmag,'all'));
ar = [range(x2),range(x3),range(x3)];

if opts.auto_lim
    lim.x = [-1,1]*max(x2)*scl.x;
    lim.y = [-1,1]*max(x3)*scl.x;
end

% ensure northbound trajectory
if not(issorted(in_situ.pos(:,2)))
    in_situ.pos = flip(in_situ.pos);
    in_situ.flow = flip(in_situ.flow);
end

% unpack in situ data + image
if strcmp(opts.pos_type,"angular")
    mlon_traj = smoothdata(in_situ.pos(:,1),"loess",16);
    mlat_traj = smoothdata(in_situ.pos(:,2),"loess",16);
    x2_traj = mlon_to_x2(mlon_traj);
    x3_traj = mlat_to_x3(mlat_traj);

    mlon_imag = image.pos(:,:,1);
    mlat_imag = image.pos(:,:,2);
    X2_imag = mlon_to_x2(mlon_imag);
    X3_imag = mlat_to_x3(mlat_imag);
else
    x2_traj = smoothdata(in_situ.pos(:,1),"loess",16);
    x3_traj = smoothdata(in_situ.pos(:,2),"loess",16);
    X2_imag = image.pos(:,:,1);
    X3_imag = image.pos(:,:,2);
end

v2_traj = in_situ.flow(:,1);
v3_traj = in_situ.flow(:,2);

% extrapolate flow data
n_ext = 32;
x3_traj_old = x3_traj;
fv2 = scatteredInterpolant(x2_traj,x3_traj,v2_traj,'natural','nearest');
fv3 = scatteredInterpolant(x2_traj,x3_traj,v3_traj,'natural','nearest');
x2_traj = [...
    x2_traj(1) + (x2_traj(2)-x2_traj(1))*(-n_ext:-1), ...
    x2_traj', ...
    x2_traj(end) + (x2_traj(end)-x2_traj(end-1))*(1:n_ext) ...
    ]';
x3_traj = [...
    x3_traj(1) + (x3_traj(2)-x3_traj(1))*(-n_ext:-1), ...
    x3_traj', ...
    x3_traj(end) + (x3_traj(end)-x3_traj(end-1))*(1:n_ext)]';
v2_traj = fv2(x2_traj,x3_traj);
v3_traj = fv3(x2_traj,x3_traj);

[~,sort_ids] = sort(x2_traj); % used for griddedInterpolant

Q = image.flux;
E0 = image.energy;
if isfield(image,'pedersen')
    fprintf('Pedersen conductance found.\n')
    SIGP = image.pedersen;
else
    Ebar = E0/1e3;
    SIGP = 40*Ebar.*sqrt(Q)./(16+Ebar.^2); % Robinson et al. (1987), Eq. (3)
end
x2_imag = X2_imag(:,1)';
x3_imag = X3_imag(1,:);

if strcmp(opts.arc_definition,'conductance')
    arc = SIGP.^2;
    ap = 1/2;
    scl.arc = 1e0;
    unt.arc = 'S';
    clm.arc = 'L18';
    lbl.arc = sprintf('Conductance (%s)',unt.arc);
elseif strcmp(opts.arc_definition,'flux')
    arc = Q;
    ap = 1;
    scl.arc = 1e0;
    unt.arc = 'mW/m^2';
    clm.arc = 'L19';
    lbl.arc = sprintf('Total energy flux (%s)',unt.arc);
end

%% calculate boundaries
num_bounds = 2;
edges = jules.tools.find_max_edges(arc,theta=0);
bsw = opts.boundary_smoothing_window;
fprintf('Boundary smoothing window is approximatly %.0f meters.\n', ...
    mean(dx3)*double(bsw))
if strcmp(opts.edge_method,'sobel')
    bounds = nan(num_bounds,size(edges,1));
    for i = 1:size(edges,1)
        [~,bounds(:,i)] = jules.tools.peak_detect(edges(i,:), ...
            num=num_bounds,smoothness=0.009);
    end
    x2_bounds_A = sort(x2_imag(2:end-1));
    x3_bounds_A = smoothdata(x3_imag(bounds(1,:)+1),2,'gaussian',bsw)';
    x2_bounds_B = x2_bounds_A;
    x3_bounds_B = smoothdata(x3_imag(bounds(2,:)+1),2,'gaussian',bsw)';
elseif strcmp(opts.edge_method,'contour')
    if all(isnan(opts.contour_values))
        edge_id = round(size(edges,1)/2);
        [~,cntr_ids] = jules.tools.peak_detect(edges(edge_id,:), ...
            num=num_bounds,smoothness=0.009);
        cntr_vals = arc(edge_id,cntr_ids);
    else
        cntr_vals = opts.contour_values.^(1/ap);
    end
    cntr_A = contour(X2_imag,X3_imag,arc,[1,1]*cntr_vals(2));
    cntr_B = contour(X2_imag,X3_imag,arc,[1,1]*cntr_vals(1));
    close(gcf)
    % lc_A = cntr_A(2,1);
    % lc_B = cntr_B(2,1);
    fprintf('Primary contour line at %.2f %s\n', ...
        cntr_A(1,1)^ap*scl.arc,unt.arc)
    fprintf('Secondary contour line at %.2f %s\n', ...
        cntr_B(1,1)^ap*scl.arc,unt.arc)
    [x2_bounds_A,x3_bounds_A] = jules.tools.get_contour(cntr_A,rank=2);
    [x2_bounds_B,x3_bounds_B] = jules.tools.get_contour(cntr_B,rank=1);
    x2_bounds_A = smoothdata(x2_bounds_A,"gaussian");
    x3_bounds_A = smoothdata(x3_bounds_A,"gaussian",bsw);
    x2_bounds_B = smoothdata(x2_bounds_B,"gaussian");
    x3_bounds_B = smoothdata(x3_bounds_B,"gaussian",bsw);
    % x2_bounds_A = smoothdata(cntr_A(1,lc_A+3:end),"gaussian");
    % x3_bounds_A = smoothdata(cntr_A(2,lc_A+3:end),"gaussian",bsw);
    % x2_bounds_A = smoothdata(cntr_A(1,2:lc_A+1),"gaussian");
    % x3_bounds_A = smoothdata(cntr_A(2,2:lc_A+1),"gaussian",bsw);
    % x2_bounds_B = smoothdata(cntr_B(1,2:lc_B+1),"gaussian");
    % x3_bounds_B = smoothdata(cntr_B(2,2:lc_B+1),"gaussian",bsw);
    [~,sort_ids_A] = sort(x2_bounds_A);
    [~,sort_ids_B] = sort(x2_bounds_B);
    x2_bounds_A = jules.tools.minsmooth(x2_bounds_A(sort_ids_A));
    x3_bounds_A = x3_bounds_A(sort_ids_A);
    x2_bounds_B = jules.tools.minsmooth(x2_bounds_B(sort_ids_B));
    x3_bounds_B = x3_bounds_B(sort_ids_B);
end

bound.A = griddedInterpolant(x2_bounds_A,x3_bounds_A);
bound.B = griddedInterpolant(x2_bounds_B,x3_bounds_B);
angle = griddedInterpolant(x2_bounds_A(2:end), ...
    atan2(diff(x3_bounds_A)',diff(x2_bounds_A)));
bound_pts = linspace(min(x2_imag),max(x2_imag),length(x2_imag));

if opts.swap_primary
    bound_A_tmp = bound.A;
    bound.A = bound.B;
    bound.B = bound_A_tmp;
end
if opts.save_data
    save('data\boundaries.mat','bound')
end

if opts.show_plots || opts.save_data
    figure
    set(gcf,'PaperPosition',[0,0,13.2,4.8])
    tiledlayout(1,2)
    ltr = opts.starting_letter;

    nexttile
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
    hold on
    pcolor(X2_imag*scl.x,X3_imag*scl.x,arc.^ap*scl.arc)
    plot(bound_pts*scl.x,bound.A(bound_pts)*scl.x,'k')
    plot(bound_pts*scl.x,bound.B(bound_pts)*scl.x,'--k')
    colormap(colorcet(clm.arc))
    clb = colorbar;
    clb.Label.String = lbl.arc;
    clb.Location = 'northoutside';
    xlabel(lbl.x); ylabel(lbl.y)
    xlim(lim.x); ylim(lim.y)
    legend('','Primary','Secondary' ...
        ,'Location','northwest','Orientation','horizontal')
    legend('boxoff')
    pbaspect(ar)

    nexttile
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8)
    hold on
    pcolor(X2_imag(2:end-1,2:end-1)*scl.x,X3_imag(2:end-1,2:end-1)*scl.x,edges)
    plot(bound_pts*scl.x,bound.A(bound_pts)*scl.x,'k')
    plot(bound_pts*scl.x,bound.B(bound_pts)*scl.x,'--k')
    colormap(colorcet(clm.arc))
    clb = colorbar;
    clb.Label.String = 'Sobel edges (a.u.)';
    clb.Location = 'northoutside';
    xlabel(lbl.x)
    yticks([])
    xlim(lim.x); ylim(lim.y)
    legend('','Primary','Secondary' ...
        ,'Location','northwest','Orientation','horizontal')
    legend('boxoff')
    pbaspect(ar)

    if all(strcmp([opts.arc_definition,opts.edge_method],["conductance","contour"]))
        tmp_fn = 'data\cond_cont.mat';
    elseif all(strcmp([opts.arc_definition,opts.edge_method],["flux","sobel"]))
        tmp_fn = 'data\flux_grad.mat';
    else
        tmp_fn = '';
    end
    if not(isempty(tmp_fn)) && opts.save_data
        fprintf('Saving %s\n',tmp_fn)
        save(tmp_fn,'X2_imag','X3_imag','edges','arc','bound','bound_pts', ...
            'scl','clm','lbl','lim','ar','ap')
    end
    if ~opts.show_plots; close all; end
end

if opts.show_plots
    figure
    hold on
    contour(X2_imag*scl.x,X3_imag*scl.x,arc.^ap*scl.arc)
    plot(bound_pts*scl.x,bound.A(bound_pts)*scl.x,'k')
    plot(bound_pts*scl.x,bound.B(bound_pts)*scl.x,'--k')
    colormap(colorcet(clm.arc))
    clb = colorbar;
    clb.Label.String = lbl.arc;
    xlabel(lbl.x); ylabel(lbl.y)
    xlim(lim.x); ylim(lim.y)
    legend('','Primary','Secondary' ...
        ,'Location','northwest','Orientation','horizontal')
    pbaspect(ar)
end

%% replicate in situ flow data
% 0 = original, 1 = replicated, a = primary boundary, b = secondary boundary
% minmooth needed for unique independent input of griddedinterpolant
% position where original trajectory meets primary boundary
traj0 = griddedInterpolant(jules.tools.minsmooth(x2_traj(sort_ids)), ...
    x3_traj(sort_ids));

x0a = fzero(@(x)(traj0(x)-bound.A(x)),0);
y0a = traj0(x0a);
% position where original trajectory meets secondary boundary
x0b = fzero(@(x)(traj0(x)-bound.B(x)),0);
y0b = traj0(x0b);
width0 = sqrt((x0b-x0a)^2+(y0b-y0a)^2);

% remove background flow
fv2_traj = griddedInterpolant(x3_traj,v2_traj);
fv3_traj = griddedInterpolant(x3_traj,v3_traj);
v2_traj_0a = fv2_traj(y0a);
v3_traj_0a = fv3_traj(y0a);
chi0 = angle(x0a);
if all(isnan(opts.flow_bg))
    alpha0 = atan2(v3_traj_0a,v2_traj_0a);
    gamma0 = pi + chi0 - alpha0;
    v2_traj_0a_rot = cos(gamma0)*v2_traj_0a - sin(gamma0)*v3_traj_0a;
    v3_traj_0a_rot = sin(gamma0)*v2_traj_0a + cos(gamma0)*v3_traj_0a;
    v_bg = [v2_traj_0a-v2_traj_0a_rot,v3_traj_0a-v3_traj_0a_rot];
    fprintf('Background flow %.2f m/s east and %.2f m/s north\n', ...
        v_bg(1),v_bg(2))
else
    v_bg = opts.flow_bg;
end

fsm = opts.flow_smoothing_window;
dx_traj = mean(sqrt(diff(x2_traj).^2+diff(x3_traj).^2));
sample_freq = 1/dx_traj;
pass_freq = sample_freq/double(fsm);
fprintf('Flow smoothing window is approximatly %.0f meters.\n',1/pass_freq)
v2_traj = smoothdata(v2_traj - v_bg(1),'gaussian',fsm);
v3_traj = smoothdata(v3_traj - v_bg(2),'gaussian',fsm);
% v2_traj = lowpass(v2_traj - v_bg(1),pass_freq,sample_freq);
% v3_traj = lowpass(v3_traj - v_bg(2),pass_freq,sample_freq);
% fv2_traj = griddedInterpolant(x3_traj,v2_traj,'spline');
% fv3_traj = griddedInterpolant(x3_traj,v3_traj,'spline');
% x2_traj = linspace(min(x2_traj),max(x2_traj),traj_us*numel(x2_traj));
% x3_traj = linspace(min(x3_traj),max(x3_traj),traj_us*numel(x3_traj));
% v2_traj = fv2_traj(x3_traj);
% v3_traj = fv3_traj(x3_traj);
% [~,sort_ids] = sort(x2_traj);

% replicate
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
    % position where replicated trajectory meets primary boundary
    traj1 = griddedInterpolant(jules.tools.minsmooth(x2_traj_tra(sort_ids)), ...
        x3_traj_tra(sort_ids));
    x1a = fzero(@(x)(traj1(x)-bound.A(x)),0);
    y1a = traj1(x1a);
    % position where replicated trajectory meets secondary boundary
    x1b = fzero(@(x)(traj1(x)-bound.B(x)),0);
    y1b = traj1(x1b);
    width1 = sqrt((x1b-x1a)^2+(y1b-y1a)^2);
    beta1 = atan2(x1b-x1a,y1b-y1a); % angle b/w line1 and vertical

    if opts.do_scale
        % rotate about (x1a,y1a) by beta
        x2_traj_rot = cos(beta1)*(x2_traj_tra - x1a) ...
            - sin(beta1)*(x3_traj_tra - y1a) + x1a;
        x3_traj_rot = sin(beta1)*(x2_traj_tra - x1a) ...
            + cos(beta1)*(x3_traj_tra - y1a) + y1a;

        % scale about p1a
        scale = width1/width0;
        x2_traj_scl = x2_traj_rot;
        x3_traj_scl = scale*(x3_traj_rot - y1a) + y1a;

        % rotate back
        x2_traj_rep(i,:) = cos(-beta1)*(x2_traj_scl - x1a) ...
            - sin(-beta1)*(x3_traj_scl - y1a) + x1a;
        x3_traj_rep(i,:) = sin(-beta1)*(x2_traj_scl - x1a) ...
            + cos(-beta1)*(x3_traj_scl - y1a) + y1a;
    else
        x2_traj_rep(i,:) = x2_traj_tra;
        x3_traj_rep(i,:) = x3_traj_tra;
    end

    if opts.do_rotate
        % rotate flows to be tangent to bound_a
        chi1 = angle(x1a) - chi0;
        v2_traj_rep(i,:) = cos(chi1)*v2_traj - sin(chi1)*v3_traj;
        v3_traj_rep(i,:) = sin(chi1)*v2_traj + cos(chi1)*v3_traj;
    else
        v2_traj_rep(i,:) = v2_traj;
        v3_traj_rep(i,:) = v3_traj;
    end

    if i == round(opts.num_replications*0.21)
        i_p1 = i;
        x2_traj_tra_p1 = x2_traj_tra; x3_traj_tra_p1 = x3_traj_tra;
        x1a_p1 = x1a; y1a_p1 = y1a;
    end
    if i == round(opts.num_replications*0.6)
        i_p2 = i;
        x2_traj_tra_p2 = x2_traj_tra; x3_traj_tra_p2 = x3_traj_tra;
        x1a_p2 = x1a; y1a_p2 = y1a;
    end
end

if opts.save_data
    save('data\reps.mat','x2_traj_rep','x3_traj_rep','v2_traj_rep','v3_traj_rep')
end

if opts.show_plots || opts.save_plots(1)
    figure
    set(gcf,'PaperPosition',[0,0,13.2,3.7]) %2.4
    tiledlayout(1,2)
    ltr = opts.starting_letter;

    [~,j_p] = min(abs(x3_traj-bound.B(x0b)));

    nexttile
    ms = 400;
    text(0.04-0.01,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
    hold on
    pcolor(X2_imag*scl.x,X3_imag*scl.x,arc.^ap*scl.arc)
    plot(bound_pts*scl.x,bound.A(bound_pts)*scl.x,'k')
    plot(bound_pts*scl.x,bound.B(bound_pts)*scl.x,'--k')
    plot(x2_traj*scl.x,x3_traj*scl.x,'r')
    plot(x2_traj_tra_p1*scl.x,x3_traj_tra_p1*scl.x,'--r')
    plot(x2_traj_tra_p2*scl.x,x3_traj_tra_p2*scl.x,'--r')
    scatter(x2_traj_tra_p1(j_p)*scl.x,x3_traj_tra_p1(j_p)*scl.x,ms,'xr')
    scatter(x2_traj_tra_p2(j_p)*scl.x,x3_traj_tra_p2(j_p)*scl.x,ms,'xr')
    plot(x2_traj_rep(i_p1,:)*scl.x,x3_traj_rep(i_p1,:)*scl.x,'b')
    plot(x2_traj_rep(i_p2,:)*scl.x,x3_traj_rep(i_p2,:)*scl.x,'b')
    scatter(x2_traj_rep(i_p1,j_p)*scl.x,x3_traj_rep(i_p1,j_p)*scl.x,ms,'xb')
    scatter(x2_traj_rep(i_p2,j_p)*scl.x,x3_traj_rep(i_p2,j_p)*scl.x,ms,'xb')
    scatter(x0a*scl.x,y0a*scl.x,ms,'xk')
    scatter(x0b*scl.x,y0b*scl.x,ms,'xr')
    scatter(x1a_p1*scl.x,y1a_p1*scl.x,ms,'xk')
    scatter(x1a_p2*scl.x,y1a_p2*scl.x,ms,'xk')
    xlim(lim.x); ylim(lim.y)
    xlabel(lbl.x); ylabel(lbl.y)
    colormap(colorcet(clm.arc))
    pbaspect(ar)

    nexttile
    vs = 10;
    cad = 16;
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8,'Color','w')
    hold on
    pcolor(X2_imag*scl.x,X3_imag*scl.x,arc.^ap*scl.arc)
    plot(bound_pts*scl.x,bound.A(bound_pts)*scl.x,'k')
    plot(bound_pts*scl.x,bound.B(bound_pts)*scl.x,'--k')
    quiver(x2_traj_rep(1:cad:end)*scl.x,x3_traj_rep(1:cad:end)*scl.x, ...
        v2_traj_rep(1:cad:end)*scl.v*vs,v3_traj_rep(1:cad:end)*scl.v*vs,0,'.-b')
    quiver(x2_traj*scl.x,x3_traj*scl.x, ...
        v2_traj*scl.v*vs,v3_traj*scl.v*vs,0,'.-r')
    % scatter(x1a*scl.x,y1a*scl.x,100,'xm')
    xlim(lim.x); ylim(lim.y)
    yticks([])
    xlabel(lbl.x)
    colormap(colorcet(clm.arc))
    pbaspect(ar)
    rectangle('Position',[73,-56,47,10],'FaceColor','w','EdgeColor','w')
    quiver(84,-51,-1*vs,0,0,'.-r','LineWidth',lw)
    text(85,-51,sprintf('%.0f %s',1,unt.v),'FontSize',fts)

    if opts.save_plots(1)
        filename = sprintf('replications%s.png',opts.suffix);
        fprintf('Saving %s\n',filename)
        saveas(gcf,fullfile(opts.direc,filename))
    end
end

if opts.show_plots
    figure
    hold on
    plot(x3_traj*scl.x,v2_traj*scl.v,'r')
    plot(x3_traj*scl.x,v3_traj*scl.v,'b')
    plot(x3_traj_old*scl.x,(in_situ.flow(:,1) - v_bg(1))*scl.v,'Color',[1,0.5,0.5])
    plot(x3_traj_old*scl.x,(in_situ.flow(:,2) - v_bg(2))*scl.v,'Color',[0.5,0.5,1])
    xlabel(lbl.y); ylabel(sprintf('Flow (%s)',unt.v))
    legend('Eastward','Northward','Eastward smoothed','Northward smoothed')
    pbaspect(ar)

    figure
    hold on
    plot(x3_traj(1:end-1)*scl.x,diff(v2_traj)*scl.v,'r')
    plot(x3_traj(1:end-1)*scl.x,diff(v3_traj)*scl.v,'b')
    xlabel(lbl.y); ylabel(sprintf('dFlow (%s)',unt.dv))
    legend('Eastward','Northward')
    pbaspect(ar)
end

%% interpolate replicated flow data
v2_int = griddata(x2_traj_rep(:),x3_traj_rep(:),v2_traj_rep(:),X2,X3,'cubic');
v3_int = griddata(x2_traj_rep(:),x3_traj_rep(:),v3_traj_rep(:),X2,X3,'cubic');
fprintf('Done interpolating replicated flow.\n')
% fv2 = griddedInterpolant(X2,X3,v2_int,'linear','linear');
% fv3 = griddedInterpolant(X2,X3,v3_int,'linear','linear');
% v2_int = fv2(X2,X3);
% v3_int = fv3(X2,X3);
E2_int =  v3_int*Bmag; % v = ExB/B^2
E3_int = -v2_int*Bmag;
E2_bg =  v_bg(2)*Bmag;
E3_bg = -v_bg(1)*Bmag;

if opts.show_plots
    figure
    vs = 2;
    hold on
    pcolor(X2_imag*scl.x,X3_imag*scl.x,arc.^ap*scl.arc)
    plot(bound_pts*scl.x,bound.A(bound_pts)*scl.x,'b')
    plot(bound_pts*scl.x,bound.B(bound_pts)*scl.x,'b')
    quiver(X2*scl.x,X3*scl.x,v2_int*scl.v*vs,v3_int*scl.v*vs,0,'.-k')
    quiver(x2_traj*scl.x,x3_traj*scl.x, ...
        v2_traj*scl.v*vs,v3_traj*scl.v*vs,0,'.-g')
    xlim(lim.x); ylim(lim.y)
    xlabel(lbl.x); ylabel(lbl.y)
    colormap(colorcet(clm.arc))
    clb = colorbar;
    clb.Label.String = lbl.arc;
    pbaspect(ar)
end

%% generate potential map
% translated from Alex Mule's python code Oct 9, 2023
% extrapolate data to avoid edge effects with fast-fourier transforms
nfill = 32;
x2_ext = [x2(1)+dx2(1)*(-nfill:-1),x2,x2(end)+dx2(end)*(1:nfill)];
x3_ext = [x3(1)+dx3(1)*(-nfill:-1),x3,x3(end)+dx3(end)*(1:nfill)];
[X2_ext,X3_ext] = ndgrid(x2_ext,x3_ext);
fE2_int = griddedInterpolant(X2,X3,E2_int,'spline','nearest');
fE3_int = griddedInterpolant(X2,X3,E3_int,'spline','nearest');
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

% find best fit harmonic function
mask_A = false(size(X2)); % select data around boundary A
mask_B = false(size(X2)); % select data around boundary B
mask_t = false(size(X2)); % select data around trajectory
for i = 1:lx2
    b = bound.A(x2(i));
    db = opts.harmonic_mask(1);
    mask_A(i,:) = (x3 >= b-db & x3 < b+db);
end
for i = 1:lx2
    b = bound.B(x2(i));
    db = opts.harmonic_mask(2);
    mask_B(i,:) = (x3 >= b-db & x3 < b+db);
end
for i = 1:lx3
    b = fzero(@(x) traj0(x)-x3(i),0);
    db = opts.harmonic_mask(3);
    mask_t(:,i) = (x2' >= b-db & x2' < b+db);
end
mask = mask_A | mask_B | mask_t;
xdata = [X2(mask),X3(mask)];
Edata = [E2_int(mask),E3_int(mask)];

if opts.fit_harmonic
    fprintf('Fitting harmonic function.\n')
    [phi,harm] = jules.tools.find_harmonic(phi0,xdata,Edata,xg);
else
    [E2_0,E3_0] = gradient(-phi0,mean(dx2),mean(dx3));
    E_0 = mean([E2_0(:);E3_0(:)])-mean([E2_int(:),E3_int(:)]);
    harm = E_0(1)*X2 + E_0(2)*X3;
    phi = phi0 + harm;
end

% adding back background electric field
if opts.add_phi_background
    phi = phi - X2*E2_bg - X3*E3_bg;
end

% remove Nyquist ringing
phi = smoothdata2(phi,"gaussian",4);

% ensure first and second differences match at north and south edges
phi(:,1) = 3*phi(:,3)-2*phi(:,4);
phi(:,2) = 2*phi(:,3)-1*phi(:,4);
phi(:,end-1) = 2*phi(:,end-2)-1*phi(:,end-3);
phi(:,end-0) = 3*phi(:,end-2)-2*phi(:,end-3);

% zero average potential
phi = phi - mean(phi,'all');

if opts.show_plots
    figure
    hold on
    surface(X2*scl.x,X3*scl.x,harm*scl.p)
    contour(X2*scl.x,X3*scl.x,double(mask),[0.5,0.5],'k')
    colormap(gca,colorcet(clm.p))
    clb = colorbar;
    clb.Label.String = sprintf('Harmonic potential (%s)',unt.p);
    xlim(lim.x); ylim(lim.y);
    xlabel(lbl.x); ylabel(lbl.y)
    pbaspect([ar(1:2),1e5])
end

%% plot final flow fields
if range(dx2) < 1e-3
    [E2,E3] = gradient(-phi',mean(dx2),mean(dx3));
else
    [E2,E3] = gradient(-phi',dx2,dx3);
end
v2 = -E3'/Bmag;
v3 =  E2'/Bmag;
divv = divergence(X3,X2,v3,v2);
divv_int = divergence(X3,X2,v3_int,v2_int);
v2_err = v2-v2_int;
v3_err = v3-v3_int;

if opts.show_plots
    figure
    vs = 2;
    hold on
    pcolor(X2_imag*scl.x,X3_imag*scl.x,arc.^ap*scl.arc)
    plot(bound_pts*scl.x,bound.A(bound_pts)*scl.x,'b')
    plot(bound_pts*scl.x,bound.B(bound_pts)*scl.x,'b')
    quiver(X2*scl.x,X3*scl.x,v2*scl.v*vs,v3*scl.v*vs,0,'.-k')
    quiver(x2_traj*scl.x,x3_traj*scl.x, ...
        v2_traj*scl.v*vs,v3_traj*scl.v*vs,0,'.-g')
    xlim(lim.x); ylim(lim.y)
    xlabel(lbl.x); ylabel(lbl.y)
    colormap(colorcet(clm.arc))
    clb = colorbar;
    clb.Label.String = lbl.arc;
    pbaspect(ar)
end

if opts.show_plots || opts.save_plots(2)
    figure
    set(gcf,'PaperPosition',[0,0,13.2,6.4])
    tiledlayout(3,3);
    ltr = opts.starting_letter;

    if opts.auto_lim
        qnt = 0.99;
        max_v = quantile(abs([v2(:)+v_bg(1);v3(:)+v_bg(2)]),qnt);
        max_dv = quantile(abs([divv(:);divv_int(:)]),qnt);
        max_p = quantile(abs(phi(:)),qnt);
        lim.v = [-1,1]*max_v*scl.v;
        lim.dv = [-1,1]*max_dv*scl.dv;
        lim.p = [-1,1]*max_p*scl.p;
    end
    
    v2_int_p = (v2_int + opts.plot_bg*v_bg(1))*scl.v;
    v3_int_p = (v3_int + opts.plot_bg*v_bg(2))*scl.v;
    v2_traj_p = (v2_traj+opts.plot_bg*v_bg(1))*scl.v;
    v3_traj_p = (v3_traj+opts.plot_bg*v_bg(2))*scl.v;
    v2_p = (v2+opts.plot_bg*v_bg(1))*scl.v;
    v3_p = (v3+opts.plot_bg*v_bg(2))*scl.v;

    % row 1
    nexttile
    title('Interpolated flow')
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
    hold on
    pcolor(X2*scl.x,X3*scl.x,v2_int_p)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj_p,v3_traj_p,'.-r')
    colormap(gca,colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = lbl.vx;
    clb.Location = 'westoutside';
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([])
    ylabel(lbl.y)
    pbaspect(ar)

    nexttile
    title('Helmholtz decomposition')
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
    hold on
    pcolor(X2*scl.x,X3*scl.x,v2_p)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj_p,v3_traj_p,'.-r')
    colormap(gca,colorcet(clm.v))
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([]); yticks([])
    pbaspect(ar)

    nexttile
    title('Difference & Potential')
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
    hold on
    pcolor(X2*scl.x,X3*scl.x,v2_err*scl.v)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj_p,v3_traj_p,'.-r')
    contour(X2*scl.x,X3*scl.x,double(mask),[0.5,0.5],'k')
    colormap(gca,colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = sprintf('\\Delta %s',lbl.vx);
    xlim(lim.x); ylim(lim.y); clim(lim.v/3)
    xticks([]); yticks([])
    pbaspect(ar)

    %row 2
    nexttile
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
    hold on
    pcolor(X2*scl.x,X3*scl.x,v3_int_p)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj_p,v3_traj_p,'.-r')
    colormap(gca,colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = lbl.vy;
    clb.Location = 'westoutside';
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([])
    ylabel(lbl.y)
    pbaspect(ar)

    nexttile
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
    hold on
    pcolor(X2*scl.x,X3*scl.x,v3_p)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj_p,v3_traj_p,'.-r')
    colormap(gca,colorcet(clm.v))
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([]); yticks([])
    pbaspect(ar)

    nexttile
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
    hold on
    pcolor(X2*scl.x,X3*scl.x,v3_err*scl.v)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj_p,v3_traj_p,'.-r')
    contour(X2*scl.x,X3*scl.x,double(mask),[0.5,0.5],'k')
    colormap(gca,colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = sprintf('\\Delta %s',lbl.vy);
    xlim(lim.x); ylim(lim.y); clim(lim.v/3)
    xticks([]); yticks([])
    pbaspect(ar)

    % row 3
    nexttile
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
    hold on
    pcolor(X2*scl.x,X3*scl.x,divv_int*scl.dv)
    colormap(gca,colorcet(clm.dv))
    clb = colorbar;
    clb.Label.String = sprintf('\\nabla\\cdot\\bfv\\rm (%s)',unt.dv);
    clb.Location = 'westoutside';
    xlim(lim.x); ylim(lim.y); clim(lim.dv)
    xlabel(lbl.x); ylabel(lbl.y)
    pbaspect(ar)

    nexttile
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
    hold on
    pcolor(X2*scl.x,X3*scl.x,divv*scl.dv)
    colormap(gca,colorcet(clm.dv))
    xlim(lim.x); ylim(lim.y); clim(lim.dv)
    yticks([])
    xlabel(lbl.x)
    pbaspect(ar)

    nexttile
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8)
    hold on
    pcolor(X2*scl.x,X3*scl.x,phi*scl.p)
    colormap(gca,colorcet(clm.p))
    clb = colorbar;
    clb.Label.String = sprintf('Potential (%s)',unt.p);
    xlim(lim.x); ylim(lim.y);
    yticks([])
    xlabel(lbl.x)
    pbaspect(ar)

    if opts.save_plots(2)
        filename = sprintf('replicated_full%s.png',opts.suffix);
        fprintf('Saving %s\n',filename)
        saveas(gcf,fullfile(opts.direc,filename))
    end
end

if opts.show_plots || opts.save_plots(3)
    figure
    set(gcf,'PaperPosition',[0,0,13.2,4.2])
    tiledlayout(1,3);
    ltr = opts.starting_letter;

    % hsv plot variables
    MLAT = 90-squeeze(xg.theta)*180/pi;
    MLON = squeeze(xg.phi)*180/pi;
    ALT = xg.alt;
    V2 = permute(repmat(v2+v_bg(1)*opts.plot_bg,[1,1,size(MLAT,1)]),[3,1,2]);
    V3 = permute(repmat(v3+v_bg(2)*opts.plot_bg,[1,1,size(MLAT,1)]),[3,1,2]);
    V2_err = permute(repmat(v2_err,[1,1,size(MLAT,1)]),[3,1,2]);
    V3_err = permute(repmat(v3_err,[1,1,size(MLAT,1)]),[3,1,2]);
    mlon_ref = mean(MLON(:));
    hsv_sat = 1e3;
    eps = 0.02;
    [hsv_map_clb,~,~,hsv_alt,hsv_alt_map] = ...
        jules.tools.hsv_params(V2,V3,MLAT,MLON,ALT,300e3,mlon_ref,hsv_sat);
    [hsv_map_clb_err,~,~,hsv_alt_err,hsv_alt_map_err] = ...
        jules.tools.hsv_params(V2_err,V3_err,MLAT,MLON,ALT,300e3,mlon_ref,hsv_sat*0.3);

    if opts.auto_lim
        qnt = 0.99;
        max_v = quantile(abs([v2(:)+v_bg(1);v3(:)+v_bg(2)]),qnt);
        max_dv = quantile(abs([divv(:);divv_int(:)]),qnt);
        max_p = quantile(abs(phi(:)),qnt);
        lim.v = [-1,1]*max_v*scl.v;
        lim.dv = [-1,1]*max_dv*scl.dv;
        lim.p = [-1,1]*max_p*scl.p;
    end

    nexttile
    title('Input flow')
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8,'Color','w'); ltr = ltr +1;
    hold on
    pcolor(X2*scl.x,X3*scl.x,hsv_alt);
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj_p,v3_traj_p,'.-r')
    contour(X2_imag*scl.x,X3_imag*scl.x,arc.^ap*scl.arc,4,'--k')%,'Color',[1,1,1]*0.4)
    colormap(gca,hsv_alt_map)
    clb = colorbar;
    colormap(clb,hsv_map_clb)
    clb.Limits = [0,1];
    clb.Ticks = [0+eps,1/4,1/2,3/4,1-eps];
    clb.TickLabels = {'W','S','E','N','W'};
    clb.Label.String = ['Sat. at ',num2str(hsv_sat*scl.v),' ',unt.v];
    clb.Location = 'southoutside';
    clim([0,1])
    xlim(lim.x); ylim(lim.y);
    xlabel(lbl.x); ylabel(lbl.y)
    pbaspect(ar)

    nexttile
    title('Input - interpolated flow')
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr +1;
    hold on
    pcolor(X2*scl.x,X3*scl.x,hsv_alt_err);
    contour(X2*scl.x,X3*scl.x,double(mask),[0.5,0.5],'--k')
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj_p,v3_traj_p,'.-r')
    colormap(gca,hsv_alt_map_err)
    clb = colorbar;
    colormap(clb,hsv_map_clb_err)
    clb.Limits = [0,1];
    clb.Ticks = [0+eps,1/4,1/2,3/4,1-eps];
    clb.TickLabels = {'W','S','E','N','W'};
    clb.Label.String = ['Sat. at ',num2str(hsv_sat*scl.v*0.3),' ',unt.v];
    clb.Location = 'southoutside';
    clim([0,1])
    xlim(lim.x); ylim(lim.y);
    yticks([])
    xlabel(lbl.x)
    pbaspect(ar)

    nexttile
    title('Input potential')
    text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8)
    hold on
    pcolor(X2*scl.x,X3*scl.x,phi*scl.p)
    colormap(gca,colorcet(clm.p))
    clb = colorbar;
    clb.Label.String = sprintf('Potential (%s)',unt.p);
    clb.Location = 'southoutside';
    xlim(lim.x); ylim(lim.y);
    yticks([])
    xlabel(lbl.x)
    pbaspect(ar)

    if opts.save_plots(3)
        filename = sprintf('replicated%s.png',opts.suffix);
        fprintf('Saving %s\n',filename)
        saveas(gcf,fullfile(opts.direc,filename))
    end

end
if ~opts.show_plots
    close all
end
end