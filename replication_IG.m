%%
direc = '//dartfs-hpc/rc/lab/L/LynchK/public_html/Gemini3D/isinglass_74';
cfg = gemini3d.read.config(direc);
direc_rep = fullfile(direc,'/ext');
cfg_rep = gemini3d.read.config(direc_rep);
xg_rep = gemini3d.grid.cartesian(cfg_rep);
load('data\replicate_data_isinglass.mat','in_situ','image')

%%
[phi,mlon,mlat,E2_bg,E3_bg] = tools.replicate(in_situ,image,xg_rep, ...
    flow_smoothing_window = 16, ...
    boundary_smoothing_window = 64, ...
    show_plots = false, ...
    save_plots = [0 0 1], ...
    direc = 'plots\paper0', ...
    suffix = '', ...
    starting_letter = 'A', ...
    add_phi_background = false, ...
    fit_harmonic = true, ...
    num_replications = 512, ... %32 or 512
    arc_definition = "conductance", ...
    edge_method = "contour", ...
    do_rotate = true, ...
    do_scale = true, ...
    harmonic_mask = [8,8,20]*1e3 ...
    );
phi_old = phi;

% mlat = mlat';
% E2_bg = 0;
% E3_bg = 0;
% close all
% 
% do_plot = false;
% lx3 = xg.lx(3);
% li = 32;
% phi_new = phi_old;
% w_max = 512;
% p = 4;
% for i=1:li
%     w = 1+(w_max-1)*(i^p-1)/(li^p-1);
%     phi_line_1 = phi(:,li+1-i);
%     phi_line_2 = phi(:,lx3-li+i);
%     smooth_phi_line_1 = smoothdata(phi_line_1,"gaussian",w);
%     smooth_phi_line_2 = smoothdata(phi_line_2,"gaussian",w);
%     if do_plot
%         figure(i)
%         hold on
% %         plot(phi_line_1,'b--')
% %         plot(smooth_phi_line_1,'b')
%         plot(phi_line_2,'r--')
%         plot(smooth_phi_line_2,'r')
%     end
%     phi_new(:,li+1-i) = smooth_phi_line_1;
%     phi_new(:,lx3-li+i) = smooth_phi_line_2;
% end
% figure(99)
% tiledlayout(1,2)
% nexttile
% surf(phi_old); shading faceted; view([45,45])
% nexttile
% surf(phi_new); shading faceted; view([45,45])

% phi = phi_new;

%%
% phi_new = smoothdata(phi_old,1,"gaussian",4);
% phi_new = smoothdata(phi_new,2,"gaussian",4);
% 
% 
% p1 = diff(phi_old,2,2);
% p2 = diff(phi_new,2,2);
% 
% hold on
% plot(p1(64,:))
% plot(p2(64,:))
% 
% phi = phi_new;
% close all
% d1phi = diff(phi,1,2);
% d2phi = diff(phi,2,2);
% figure
% surface(d1phi)
% figure
% surface(10*d2phi)

%%
[X2_rep,X3_rep] = ndgrid(xg_rep.x2(3:end-2),xg_rep.x3(3:end-2));
x3_rep = xg_rep.x3(3:end-2);
clear('xg_rep')

clear('xg')
xg = gemini3d.grid.cartesian(cfg);
[X2,X3] = ndgrid(xg.x2(3:end-2),xg.x3(3:end-2));

fphi = griddedInterpolant(X2_rep,X3_rep,phi,'spline');
phi = fphi(X2,X3);

%%
close all
hold on
plot(x3_rep,phi_old(64,:))
plot(xg.x3(3:end-2),phi(64,:))

%%
% E2_bg = 0;
% E3_bg = 0;
save(fullfile(direc,'ext','potential_map.mat'),'phi','E2_bg','E3_bg')

%%
dat = gemini3d.read.frame(direc,'time',cfg.times(2));

%%
ll=5;
[X2,X3] = ndgrid(xg.x2(3+ll:end-2-ll),xg.x3(3+ll:end-2-ll));
[dX2,dX3] = ndgrid(xg.dx2h(ll+1:end-ll),xg.dx3h(ll+1:end-ll));
dA = dX2.*dX3;
FAC = -squeeze(dat.J1(end,ll+1:end-ll,ll+1:end-ll))*1e6;
sum(FAC.*dA,'all')
mean(FAC(:))
median(FAC(:))
pcolor(X2,X3,FAC); colorbar; clim([-1,1]*quantile(FAC(:),0.99)); colormap(gca,colorcet('D1A'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% *** old ***
%% init
sc1 = 10;
sc2 = 1.5;

do_plot = true;
save_plot = true;
auto_lim = true;
do_A = true;
do_B = true;
do_C = true;
%#ok<*UNRCH>

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
tools.setall(0,'FontName',ftn)
tools.setall(0,'FontSize',fts)
tools.setall(0,'Multiplier',1)

scl.x = 1e-3; scl.v = 1e-3; scl.dv = 1e0; scl.p = 1e-3;
unt.x = 'km'; unt.v = 'km/s'; unt.dv = 'Hz'; unt.p = 'kV'; unt.q = 'mW/m^2';
clm.v = 'D2'; clm.dv = 'CBD1'; clm.p = 'D10'; clm.q = 'L19';
lim.x = [-1,1]*125; lim.y = [-1,1]*59;  lim.v = [-1,1]*1.3; lim.dv = [-1,1]*0.1;

lbl.x = sprintf('M. east (%s)',unt.x);
lbl.y = sprintf('M. north (%s)',unt.x);
lbl.vx = sprintf('v_{east} (%s)',unt.v);
lbl.vy = sprintf('v_{north} (%s)',unt.v);

%% unpack grid and configuration data
direc = '//dartfs-hpc/rc/lab/L/LynchK/public_html/Gemini3D/isinglass_05';
direc = fullfile(direc);
cfg = gemini3d.read.config(direc);
xg = gemini3d.grid.cartesian(cfg);
MLAT = 90-squeeze(xg.theta(end,:,:))*180/pi;
MLON = squeeze(xg.phi(end,:,:))*180/pi;
x2 = xg.x2(3:end-2);
x3 = xg.x3(3:end-2);
dx2 = xg.dx2h;
dx3 = xg.dx3h;
lx2 = xg.lx(2); lx3 = xg.lx(3);
[X2,X3] = ndgrid(x2,x3);
[DX2,DX3] = ndgrid(dx2,dx3);
Bmag = abs(mean(xg.Bmag,'all'));
mlon_to_x2 = griddedInterpolant(MLON(:,1),xg.x2(3:end-2));
mlat_to_x3 = griddedInterpolant(MLAT(1,:),xg.x3(3:end-2));
ar = [range(x2),range(x3),range(x3)];

%% determine arc boundaries
it_part = 4; % datetime = 2017 3 2 7 52 40 + it_part-1
load(fullfile(direc,'ext','particles.mat'))
datetime_part = datetime(outputdate(it_part,:));
Q_map_part = Qit(:,:,it_part);
E0_map_part = E0it(:,:,it_part);
x2_part = mlon_to_x2(mlon);
x3_part = mlat_to_x3(mlat);
[X2_part,X3_part] = ndgrid(x2_part,x3_part);

edges = tools.find_max_edges(Q_map_part);
num_bounds = 2;
bounds = nan(num_bounds,size(edges,1));
for i = 1:size(edges,1)
    [~,bounds(:,i)] = tools.peak_detect(edges(i,:), ...
        num=num_bounds,smoothness=0.01);
end

bounds_x2 = x2_part(2:end-1);
bounds_x3 = smoothdata(x3_part(bounds+1),2,"gaussian");
bound_a = griddedInterpolant(bounds_x2,bounds_x3(1,:));
bound_b = griddedInterpolant(bounds_x2,bounds_x3(2,:));
angle = griddedInterpolant(bounds_x2(1:end-1), ...
    atan2(diff(bounds_x3(1,:)),diff(bounds_x2)));

if do_plot
    close all
    figure
    pts_bound = -155e3:5e3:155e3;
    hold on
    pcolor(X2_part*scl.x, X3_part*scl.x, Q_map_part)
    scatter(bounds_x2*scl.x, bounds_x3(1,:)*scl.x,'k')
    scatter(bounds_x2*scl.x, bounds_x3(2,:)*scl.x,'k')
    plot(pts_bound*scl.x, bound_a(pts_bound)*scl.x,'k')
    plot(pts_bound*scl.x, bound_b(pts_bound)*scl.x,'k')
    shading flat
    clb = colorbar;
    clb.Label.String = ['Total precipitating energy flux (',unt.q,')'];
    colormap(gca,colorcet(clm.q))
    xlim(lim.x)
    ylim(lim.y)
    xlabel(lbl.x)
    ylabel(lbl.y)
    pbaspect(ar)
end

clear('glon','glat','mlon','mlat','Qit','E0it','outputdate')

%% replicate rocket data
do_scale = true;
do_rotate = true;
it_isin = 15;
v2_bg = 0;
v3_bg = -500; % chi-by-eye constant background flow removal

load(fullfile(direc,'ext','flow_data.mat'))
x2_traj = mlon_to_x2(mlon);
x3_traj = mlat_to_x3(mlat);
v2_traj = smoothdata(v_geom_east - v2_bg,'gaussian',1);
v3_traj = smoothdata(v_geom_north - v3_bg,'gaussian',1);
datetime_isin = datetimes(it_isin);

fprintf('Particle data datetime: %s\n', datetime_part)
fprintf('Rocket data datetime: %s\n', datetime_isin)

% 0 = original, 1 = replicated, a = primary bound, b = secondary bound
traj0 = griddedInterpolant(flip(x2_traj),flip(x3_traj));
x0a = fzero(@(x)(traj0(x)-bound_a(x)),0);
x0b = fzero(@(x)(traj0(x)-bound_b(x)),0);
y0a = traj0(x0a);
y0b = traj0(x0b);
width0 = sqrt((x0b-x0a)^2+(y0b-y0a)^2);
beta = atan2(x0b-x0a,y0b-y0a); % angle b/w bound-traj inters. and vertical

num_reps = 256;
dxs = linspace(-120e3,190e3,num_reps); % eastward displacements
x2_traj_rep = nan(length(dxs),length(x2_traj));
x3_traj_rep = nan(length(dxs),length(x3_traj));
v2_traj_rep = nan(length(dxs),length(v2_traj));
v3_traj_rep = nan(length(dxs),length(v3_traj));
for i = 1:length(dxs)
    % translate
    dx = dxs(i);
    dy = bound_a(x0a+dx)-y0a;

    x2_traj_tra = x2_traj + dx;
    x3_traj_tra = x3_traj + dy;

    % determine width at position 1
    traj1 = griddedInterpolant(flip(x2_traj_tra),flip(x3_traj_tra));
    x1a = fzero(@(x)(traj1(x)-bound_a(x)),0);
    x1b = fzero(@(x)(traj1(x)-bound_b(x)),0);
    y1a = traj1(x1a);
    y1b = traj1(x1b);
    width1 = sqrt((x1b-x1a)^2+(y1b-y1a)^2);

    if do_scale
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
    
    if do_rotate
        % rotate flows to be tangent to bound_a
        alpha = angle(x1a);
        v2_traj_rep(i,:) = cos(alpha)*v2_traj - sin(alpha)*v3_traj;
        v3_traj_rep(i,:) = sin(alpha)*v2_traj + cos(alpha)*v3_traj;
    else
        v2_traj_rep(i,:) = v2_traj;
        v3_traj_rep(i,:) = v3_traj;
    end

    if i==50 && do_plot
        close all
        figure
        hold on
        pcolor(X2_part*scl.x,X3_part*scl.x,Q_map_part); shading flat
        colormap(gca,"gray")
        plot(bounds_x2*scl.x,bounds_x3(1,:)*scl.x,'b')
        plot(bounds_x2*scl.x,bounds_x3(2,:)*scl.x,'b')
        scatter(x0a*scl.x,y0a*scl.x,50,'Filled','y')
        scatter(x0b*scl.x,y0b*scl.x,50,'Filled','y')
        scatter(x1a*scl.x,y1a*scl.x,50,'Filled','y')
        scatter(x1b*scl.x,y1b*scl.x,50,'Filled','y')
        quiver(x2_traj*scl.x,x3_traj*scl.x, ...
            v2_traj*scl.v*sc1,v3_traj*scl.v*sc1,0,'.-g')
        plot(x2_traj_tra*scl.x,x3_traj_tra*scl.x,'b')
        plot(x2_traj_rot*scl.x,x3_traj_rot*scl.x,'r')
        plot(x2_traj_scl*scl.x,x3_traj_scl*scl.x,'b--')
        quiver(x2_traj_rep(i,:)*scl.x,x3_traj_rep(i,:)*scl.x, ...
            v2_traj*scl.v*sc1,v3_traj*scl.v*sc1,0,'.-r')
        xlim(lim.x)
        ylim(lim.y)
        xlabel(lbl.x)
        ylabel(lbl.y)
        pbaspect(ar)
    end
end

if do_plot
    figure
    hold on
    pcolor(X2_part*scl.x,X3_part*scl.x,Q_map_part); shading flat
    colormap(gca,"gray")
    plot(bounds_x2*scl.x,bounds_x3(1,:)*scl.x,'b')
    plot(bounds_x2*scl.x,bounds_x3(2,:)*scl.x,'b')
    quiver(x2_traj_rep*scl.x,x3_traj_rep*scl.x, ...
        v2_traj_rep*scl.v,v3_traj_rep*scl.v,0,'.-r')
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,0,'.-g')
    xlim(lim.x)
    ylim(lim.y)
    xlabel(lbl.x)
    ylabel(lbl.y)
    pbaspect(ar)

    figure
    hold on
    plot(time,v_geom_east - v2_bg,'Color',[1,0.4,0.4])
    plot(time,v_geom_north - v3_bg,'Color',[0.4,0.4,1])
    plot(time,v2_traj,'r')
    plot(time,v3_traj,'b')
    xlabel('Flight time (s)')
    ylabel('Flow (m/s)')
    legend('east','east smoothed','north','north smoothed')
    pbaspect(ar)
end

clear('lon_100km','lat_100km','mlon','mlat', ...
    'v_geom_east','v_geom_north','v_geog_north','v_geog_north', ...
    'time','datetimes','comment','note')

%% interpolate replicated flow data
fv2 = scatteredInterpolant(x2_traj_rep(:),x3_traj_rep(:),v2_traj_rep(:));
fv3 = scatteredInterpolant(x2_traj_rep(:),x3_traj_rep(:),v3_traj_rep(:));
v2_int = fv2(X2,X3);
v3_int = fv3(X2,X3);
E2_int =  v3_int*Bmag; % v = ExB/B^2
E3_int = -v2_int*Bmag;
E2_bg =  v3_bg*Bmag;
E3_bg = -v2_bg*Bmag;

if do_plot
    close all
    figure
    hold on
    pcolor(X2*scl.x,X3*scl.x,v2_int*scl.v); shading flat
    colorbar; colormap(gca,colorcet(clm.v)); clim(lim.v)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v*sc1,v3_traj*scl.v*sc1,0,'.-k')
    xlim(lim.x)
    ylim(lim.y)
    xlabel(lbl.x)
    ylabel(lbl.y)
    pbaspect(ar)

    figure
    hold on
    pcolor(X2*scl.x,X3*scl.x,v3_int*scl.v); shading flat
    colorbar; colormap(gca,colorcet(clm.v)); clim(lim.v)
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v*sc1,v3_traj*scl.v*sc1,0,'.-k')
    xlim(lim.x)
    ylim(lim.y)
    xlabel(lbl.x)
    ylabel(lbl.y)
    pbaspect(ar)

    figure
    vs = 2;
    hold on
    pcolor(X2_part*scl.x,X3_part*scl.x,Q_map_part); shading flat
    colormap(gca,colorcet(clm.q))
    plot(bounds_x2*scl.x,bounds_x3(1,:)*scl.x,'b')
    plot(bounds_x2*scl.x,bounds_x3(2,:)*scl.x,'b')
    quiver(X2*scl.x,X3*scl.x,v2_int*scl.v*vs,v3_int*scl.v*vs,0,'.-r')
    quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v*vs,v3_traj*scl.v*vs,0,'.-g')
    xlim(lim.x)
    ylim(lim.y)
    xlabel(lbl.x)
    ylabel(lbl.y)
    pbaspect(ar)
end

%% generate potential map - option A
if do_A
    pseudo_basis = false;
    usepar = false;
    look_for_file = true;
    suf = 'v7';
    dec = [1,1];
    numf = 64;
    
    if pseudo_basis
        bag_x2 = X2(1:dec(1):end,1:dec(2):end);
        bag_x3 = X3(1:dec(1):end,1:dec(2):end);
        bag_v2 = v2_int(1:dec(1):end,1:dec(2):end);
        bag_v3 = v3_int(1:dec(1):end,1:dec(2):end);
        bag_x = [bag_x2(:),bag_x3(:)];
        bag_v = [bag_v2(:),bag_v3(:)];
        boundary = [bounds_x2;bounds_x3(1,:)]';
    
        if do_plot
            close all
    
            figure
            hold on
            pcolor(X2*scl.x,X3*scl.x,v2_int*scl.v); shading flat
            colorbar; colormap(gca,clm.v); clim(lim.v)
            quiver(bag_x2*scl.x,bag_x3*scl.x,bag_v2*scl.v,bag_v3*scl.v,0,'.-k')
            pbaspect(ar)
            cont = input('Looks good? (y/n) ',"s");
        else
            cont = 'y'
        end
    
        mat_fn = sprintf('phitopA_%i_%i_%i_%s.mat',dec(1),dec(2),numf,suf);
        mat_fn = fullfile('data',mat_fn);
        if strcmp(cont,'y')
            if exist(mat_fn,'file') == 2 && look_for_file
                fprintf('File found: %s\n',mat_fn)
                load(mat_fn)
            else
                phitopA = tools.reconstruct(bag_x,bag_v,boundary,xg,numf=numf, ...
                    usepar=usepar);
                save(mat_fn,'phitopA')
            end
        end
    else
        mat_fn = sprintf('phitopA_%i_%i_%s.mat',dec(1),dec(2),suf);
        mat_fn = fullfile('data',mat_fn);
        if exist(mat_fn,'file') == 2 && look_for_file
            fprintf('File found: %s\n',mat_fn)
            load(mat_fn)
        else
    %         phi_0 = zeros(size(E2));
    %         for ix2 = 1:lx2
    %             for ix3 = 1:lx3
    %                 phi_0(ix2,ix3) = -( ...
    %                     dot(E2(1:ix2,1), dx2(1:ix2)) + ...
    %                     dot(E3(ix2,1:ix3), dx3(1:ix3)) ...
    %                     );
    %             end
    %         end
    %         phi_0 = phi_0 - mean(phi_0,'all');
    
            phitopA = tools.basic_reconstruct(v2_int,v3_int,xg ...
                ,usepar=usepar,decimation=dec,ic_iteration=true);
            save(mat_fn,'phitopA','dec','usepar')
        end
    end
    phitopA = phitopA - mean(phitopA,'all');
end
%% generate potential map - option B
if do_B
    look_for_file = true;
    suf = 'v1';
    ns = [4,256];
    
    [s2s,s3s] = ndgrid(round(linspace(1,lx2,ns(1))),round(linspace(1,lx3,ns(2))));
    mat_fn = sprintf('phitopB_%i_%i_%s.mat',ns(1),ns(2),suf);
    mat_fn = fullfile('data',mat_fn);
    if exist(mat_fn,'file') == 2 && look_for_file
        fprintf('File found: %s\n',mat_fn)
        load(mat_fn)
    else
        phitopB = zeros([size(E2_int),prod(ns)]);
        for n = 1:ns(1)*ns(2)
            s2 = s2s(n);
            s3 = s3s(n);
            for ix2 = 1:lx2
                if ix2 < s2
                    d2 = -1;
                else
                    d2 = 1;
                end
                for ix3 = 1:lx3
                    if ix3 < s3
                        d3 = -1;
                    else
                        d3 = 1;
                    end
                    phitopB(ix2,ix3,n) = -( ...
                        dot(E2_int(s2:d2:ix2,s3), dx2(s2:d2:ix2)*d2) ...
                        + dot(E3_int(ix2,s3:d3:ix3), dx3(s3:d3:ix3)*d3) ...
                        );
                end
            end
            %     figure(i)
            %     hold on
            %     pcolor(X2,X3,phitopB); shading flat
        end
        phitopB = mean(phitopB,3);
        phitopB = phitopB - mean(phitopB,'all');
        save(mat_fn,'phitopB','ns')
    end
end

%% generate potential map - option C
% translated from Alex Mule's python code Oct 9, 2023

if do_C
    % extrapolate data to avoid edge effects with fourier transforms
    nfill = 16;
    x2_ext = [x2(1)+dx2(1)*(-nfill:-1),x2,x2(end)+dx2(end)*(1:nfill)];
    x3ext = [x3(1)+dx3(1)*(-nfill:-1),x3,x3(end)+dx3(end)*(1:nfill)];
    [X2_ext,X3_ext] = ndgrid(x2_ext,x3ext);
    fE2_map = griddedInterpolant(X2,X3,E2_int,'linear','nearest');
    fE3_map = griddedInterpolant(X2,X3,E3_int,'linear','nearest');
    E2_map_ext = fE2_map(X2_ext,X3_ext);
    E3_map_ext = fE3_map(X2_ext,X3_ext);
    
    % wave vector convention: 2 pi ( -f_Ny : +f_Ny )
    k2_Ny = 2*pi/(2*mean(dx2));
    k3_Ny = 2*pi/(2*mean(dx3));
    k2 = linspace(-k2_Ny,k2_Ny,lx2+2*nfill);
    k3 = linspace(-k3_Ny,k3_Ny,lx3+2*nfill);
    [K2,K3] = ndgrid(k2,k3);
    Khat2 = K2./sqrt(K2.^2 + K3.^2);
    Khat3 = K3./sqrt(K2.^2 + K3.^2);
    
    % testing matlab fft2 and wave-vector convention
    % k2_target = 4.2*2*pi/5e4;
    % k3_target = 5.4*2*pi/5e4;
    % test_image = cos(X2*k2_target) + cos(X3*k3_target);
    % test_fft = fftshift(fft2(test_image));
    % 
    % [~,i2] = max(abs(test_fft(:,end/2+1)));
    % [~,i3] = max(abs(test_fft(end/2+1,:)));
    % k2_fft = abs(k2(i2));
    % k3_fft = abs(k3(i3));
    % error2 = 100*(1-k2_fft/k2_target);
    % error3 = 100*(1-k3_fft/k3_target);
    % fprintf('%% error in finding k2_target = %.5f\n',error2)
    % fprintf('%% error in finding k3_target = %.5f\n',error3)
    % 
    % if plotting
    %     figure
    %     hold on
    %     pcolor(X2-mean(diff(x2)),X3-mean(diff(x3)),test_image); shading flat
    %     quiver(0,0,2*pi/k2_target,0,'r')
    %     quiver(0,0,0,2*pi/k3_target,'r')
    %     
    %     figure
    %     hold on
    %     pcolor(K2-mean(diff(k2)),K3-mean(diff(k3)),abs(test_fft)); shading flat
    %     xlim([-1,1]*1e-3)
    %     ylim([-1,1]*1e-3)
    %     scatter(k2_target,0,1e3,'xr')
    %     scatter(k2_fft,0,1e3,'xg')
    %     scatter(0,k3_target,1e3,'xr')
    %     scatter(0,k3_fft,1e3,'xg')
    % end
    
    % Fourier transform of electric field
    G2 = fftshift(fft2(E2_map_ext));
    G3 = fftshift(fft2(E3_map_ext));
    
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
    phitopC = phi0 - mean(E2_int-E20,'all').*X2 - mean(E3_int-E30,'all').*X3;
    phitopC = phitopC - mean(phitopC,'all');
end

%% plot potentials A, B, and C
if do_plot
    close all
    figure
    tiledlayout(2,2)
    nexttile
    if do_A
        pcolor(X2,X3,phitopA); shading flat; colorbar; clim([-1,1]*2e3)
        pbaspect(ar)
    end
    nexttile
    if do_B
        pcolor(X2,X3,phitopB); shading flat; colorbar; clim([-1,1]*2e3)
        pbaspect(ar)
    end
    nexttile
    if do_C
        pcolor(X2,X3,phitopC); shading flat; colorbar; clim([-1,1]*2e3)
        pbaspect(ar)
        nexttile
        delta_phi = phi0-phitopA;
        surface(X2/1e3,X3/1e3,delta_phi/1e3)
        xlabel('East (km)'); ylabel('North (km)'); zlabel('Et tu brutus? - Alex (kV)')
        view(45,45)
        grid
    end
end

%% plot interpolated flow against options A, B, and C
close all

add_bg = 1;
suff = 'v1';

if ~do_A; phitopA = zeros(size(X2)); end
if ~do_B; phitopB = zeros(size(X2)); end
if ~do_C; phitopC = zeros(size(X2)); end

if range(dx2) < 1e-3
    [E2_int_A,E3_int_A] = gradient(-phitopA',mean(dx2),mean(dx3));
    [E2_int_B,E3_int_B] = gradient(-phitopB',mean(dx2),mean(dx3));
    [E2_int_C,E3_int_C] = gradient(-phitopC',mean(dx2),mean(dx3));
else
    [E2_int_A,E3_int_A] = gradient(-phitopA',dx2,dx3);
    [E2_int_B,E3_int_B] = gradient(-phitopB',dx2,dx3);
    [E2_int_C,E3_int_C] = gradient(-phitopC',dx2,dx3);
end
v2_A = -E3_int_A'/Bmag;
v2_B = -E3_int_B'/Bmag;
v2_C = -E3_int_C'/Bmag;
v3_A =  E2_int_A'/Bmag;
v3_B =  E2_int_B'/Bmag;
v3_C =  E2_int_C'/Bmag;

if do_plot
    for c = 'ABC'
        if c=='A'; v2_tmp = v2_A; v3_tmp = v3_A; end
        if c=='B'; v2_tmp = v2_B; v3_tmp = v3_B; end
        if c=='C'; v2_tmp = v2_C; v3_tmp = v3_C; end
        figure
        vs = 2;
        hold on
        pcolor(X2_part*scl.x,X3_part*scl.x,Q_map_part); shading flat
        colormap(gca,colorcet(clm.q))
        plot(bounds_x2*scl.x,bounds_x3(1,:)*scl.x,'b')
        plot(bounds_x2*scl.x,bounds_x3(2,:)*scl.x,'b')
        quiver(X2*scl.x,X3*scl.x,v2_tmp*scl.v*vs,v3_tmp*scl.v*vs,0,'.-r')
        quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v*vs,v3_traj*scl.v*vs,0,'.-g')
        xlim(lim.x)
        ylim(lim.y)
        xlabel(lbl.x)
        ylabel(lbl.y)
        pbaspect(ar)
    end
end

%%
divv_int = divergence(X3,X2,v3_int,v2_int);

if auto_lim
    qnt = 0.95;
    max_v = quantile(abs([v2_int(:)+v2_bg;v3_int(:)+v3_bg]),qnt);
    max_dv = quantile(abs(divv_int(:)),qnt);
    max_p = quantile(abs(phitopA(:)),qnt);
    lim.v = [-1,1]*max_v*scl.v;
    lim.dv = [-1,1]*max_dv*scl.dv;
    lim.p = [-1,1]*max_p*scl.p;
end

figure
set(gcf,'PaperPosition',[0,0,13.2,8]) %4.3
tiledlayout(4,3);
ltr = 65;

% row 1
nexttile
title('Interpolated flow')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2*scl.x,X3*scl.x,(v2_int+v2_bg)*scl.v)
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xticks([])
ylabel(lbl.y)
pbaspect(ar)

nexttile
title('Brute force')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2*scl.x,X3*scl.x,(v2_A+v2_bg)*scl.v)
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xticks([]); yticks([])
pbaspect(ar)

nexttile
title('Helmholtz decomposition')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2*scl.x,X3*scl.x,(v2_C+v2_bg)*scl.v)
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
clb = colorbar;
clb.Label.String = lbl.vx;
xlim(lim.x); ylim(lim.y); clim(lim.v)
xticks([]); yticks([])
pbaspect(ar)

% row 2
nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_int+v3_bg)*scl.v)
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xticks([])
ylabel(lbl.y)
pbaspect(ar)

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_A+v3_bg)*scl.v)
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xticks([]); yticks([])
pbaspect(ar)

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_C+v3_bg)*scl.v)
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
clb = colorbar;
clb.Label.String = lbl.vy;
xlim(lim.x); ylim(lim.y); clim(lim.v)
xticks([]); yticks([])
pbaspect(ar)

% row 3
nexttile
% title('Divergence & Potential')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2*scl.x,X3*scl.x,divv_int*scl.dv)
colormap(gca,colorcet(clm.dv))
clb = colorbar;
clb.Label.String = sprintf('\\nabla\\cdot\\bfv\\rm (%s)',unt.dv);
clb.Location = 'westoutside';
xlim(lim.x); ylim(lim.y); clim(lim.dv)
xticks([])
ylabel(lbl.y)
pbaspect(ar)

nexttile
% title('Brute force - Interp.')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2*scl.x,X3*scl.x,(v2_A-v2_int)*scl.v)
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v/3)
xticks([]); yticks([])
pbaspect(ar)

nexttile
% title('Helmholts - Brute force')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2*scl.x,X3*scl.x,(v2_C - v2_A)*scl.v)
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
clb = colorbar;
clb.Label.String = sprintf('\\Delta %s',lbl.vx);
xlim(lim.x); ylim(lim.y); clim(lim.v/3)
xticks([]); yticks([])
pbaspect(ar)

% row 4
nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2*scl.x,X3*scl.x,100*(phitopA-phitopC)*scl.p/lim.p(2))
colormap(gca,colorcet(clm.p))
clb = colorbar;
clb.Label.String = '\Delta\phi/\phi_{max} (%)';
clb.Location = 'westoutside';
xlim(lim.x); ylim(lim.y); %clim([-1,1]*10)
xlabel(lbl.x); ylabel(lbl.y)
pbaspect(ar)

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_A-v3_int)*scl.v)
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v/3)
yticks([])
xlabel(lbl.x);
pbaspect(ar)

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8)
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_C - v3_A)*scl.v)
quiver(x2_traj*scl.x,x3_traj*scl.x,v2_traj*scl.v,v3_traj*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
clb = colorbar;
clb.Label.String = sprintf('\\Delta %s',lbl.vy);
xlim(lim.x); ylim(lim.y); clim(lim.v/3)
yticks([])
xlabel(lbl.x);
pbaspect(ar)

% if do_A
%     if pseudo_basis
%         filename = sprintf('phitopA_%i_%i_%i_%s.png',dec(1),dec(2),numf,suf);
%     else
%         filename = sprintf('phitopA_%i_%i_%s.png',dec(1),dec(2),suf);
%     end
%     filename = fullfile('plots',filename);
% else
%     filename = fullfile('plots','phitop.png');
% end
% filename = sprintf('potential_options_%s.png', ...
    % datetime(datetime,'Format','MMdd''_''HHmmss'));
% filename = fullfile('plots',filename);
filename = fullfile('plots','paper0','potential_options.png');
if save_plot
    fprintf('Saving file:c %s\n',filename)
    saveas(gcf,filename)
end
close all

%%
close all
if auto_lim
    qnt = 0.95;
    max_v = quantile(abs([v2_int(:)+v2_bg;v3_int(:)+v3_bg]),qnt);
    max_dv = quantile(abs(divv_int(:)),qnt);
    max_p = quantile(abs(phitopA(:)),qnt);
    lim.v = [-1,1]*max_v*scl.v;
    lim.dv = [-1,1]*max_dv*scl.dv;
    lim.p = [-1,1]*max_p*scl.p;
end

% find best fit harmonic function
xdata = [X2(:),X3(:)];
Edata = [E2_int(:),E3_int(:)];
[~,harm1] = tools.find_harmonic(phi0,xdata,Edata,xg);

harmonic_mask = [1,1,2]*10e3;
mask_A = false(size(X2)); % select data around boundary A
mask_B = false(size(X2)); % select data around boundary B
mask_t = false(size(X2)); % select data around trajectory
for i = 1:lx2
    b = bound_a(x2(i));
    db = harmonic_mask(1);
    mask_A(i,:) = (x3 >= b-db & x3 < b+db);
end
for i = 1:lx2
    b = bound_b(x2(i));
    db = harmonic_mask(2);
    mask_B(i,:) = (x3 >= b-db & x3 < b+db);
end
for i = 1:lx3
    b = fzero(@(x) traj0(x)-x3(i),0);
    db = harmonic_mask(3);
    mask_t(:,i) = (x2' >= b-db & x2' < b+db);
end
mask = mask_A | mask_B | mask_t;
xdata = [X2(mask),X3(mask)];
Edata = [E2_int(mask),E3_int(mask)];
[~,harm2] = tools.find_harmonic(phi0,xdata,Edata,xg);

%%
close all

figure
set(gcf,'PaperPosition',[0,0,13.2,4.2]) %4.3
tiledlayout(1,3);
ltr = 'A';

nexttile
lc = 6;
vw = [130,20];
title('Brute force - \phi_0(x,y)')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
surface(X2(1:lc:end,1:lc:end)*scl.x,X3(1:lc:end,1:lc:end)*scl.x, ...
    (phitopA(1:lc:end,1:lc:end)-phi0(1:lc:end,1:lc:end))*scl.p,'EdgeColor','k')
colormap(gca,colorcet(clm.p))
xlim(lim.x); ylim(lim.y); zlim(1.4*lim.p); clim(lim.p)
ylbl = ylabel(lbl.y);
ylbl.Position = [180,100,-2];
zlabel(sprintf('Potential (%s)',unt.p))
pbaspect(ar)
view(vw)
grid

nexttile
title('Harmonic fit')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
surface(X2(1:lc:end,1:lc:end)*scl.x,X3(1:lc:end,1:lc:end)*scl.x, ...
    harm1(1:lc:end,1:lc:end)*scl.p,'EdgeColor','k')
colormap(gca,colorcet(clm.p))
xlim(lim.x); ylim(lim.y); zlim(1.4*lim.p); clim(lim.p)
set(gca,'ZTickLabel',[])
pbaspect(ar)
view(vw)
grid

nexttile
harm2_new = harm2;
harm2_new(not(mask)) = nan;
title('Harmonic fit (masked)')
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8)
hold on
surface(X2(1:lc:end,1:lc:end)*scl.x,X3(1:lc:end,1:lc:end)*scl.x, ...
    harm2(1:lc:end,1:lc:end)*scl.p+2,'EdgeColor','k','FaceAlpha',0.3,'EdgeAlpha',0.3)
surface(X2(1:lc:end,1:lc:end)*scl.x,X3(1:lc:end,1:lc:end)*scl.x, ...
    harm2_new(1:lc:end,1:lc:end)*scl.p+2,'EdgeColor','k')
contour(X2*scl.x,X3*scl.x,double(mask),'r')
colormap(gca,colorcet(clm.p))
xlim(lim.x); ylim(lim.y); zlim(1.4*lim.p+2); clim(lim.p+2)
xlbl = xlabel(lbl.x);
xlbl.Position = [-30,30,-2];
set(gca,'ZTickLabel',[])
pbaspect(ar)
view(vw)
grid

filename = fullfile('plots','paper0','harmonic.png');
if save_plot
    fprintf('Saving file:c %s\n',filename)
    saveas(gcf,filename)
end
close all

%%
[E2_int_A,E3_int_A] = gradient(-phitopA',mean(dx2),mean(dx3));
[E2_int_1,E3_int_1] = gradient(-(phi0+harm1)',mean(dx2),mean(dx3));
[E2_int_2,E3_int_2] = gradient(-(phi0+harm2)',mean(dx2),mean(dx3));

v2_A = -E3_int_A'/Bmag;
v2_1 = -E3_int_1'/Bmag;
v2_2 = -E3_int_2'/Bmag;
v3_A =  E2_int_A'/Bmag;
v3_1 =  E2_int_1'/Bmag;
v3_2 =  E2_int_2'/Bmag;

fv2_int_A = griddedInterpolant(X2,X3,v2_A);
fv2_int_1 = griddedInterpolant(X2,X3,v2_1);
fv2_int_2 = griddedInterpolant(X2,X3,v2_2);
fv3_int_A = griddedInterpolant(X2,X3,v3_A);
fv3_int_1 = griddedInterpolant(X2,X3,v3_1);
fv3_int_2 = griddedInterpolant(X2,X3,v3_2);

v2_bnd_A = fv2_int_A(bounds_x2,bounds_x3(1,:));
v2_bnd_1 = fv2_int_1(bounds_x2,bounds_x3(1,:));
v2_bnd_2 = fv2_int_2(bounds_x2,bounds_x3(1,:));
v3_bnd_A = fv3_int_A(bounds_x2,bounds_x3(1,:));
v3_bnd_1 = fv3_int_1(bounds_x2,bounds_x3(1,:));
v3_bnd_2 = fv3_int_2(bounds_x2,bounds_x3(1,:));

dangle_A = rad2deg(atan(v3_bnd_A./v2_bnd_A)-angle(bounds_x2));
dangle_1 = rad2deg(atan(v3_bnd_1./v2_bnd_1)-angle(bounds_x2));
dangle_2 = rad2deg(atan(v3_bnd_2./v2_bnd_2)-angle(bounds_x2));

disp('----------------')
disp(quantile(dangle_A,[0.1,0.5,0.9]))
disp(quantile(dangle_1,[0.1,0.5,0.9]))
disp(quantile(dangle_2,[0.1,0.5,0.9]))

close all

figure
hold on
histogram(dangle_A,32,'FaceColor','k','Normalization','probability')
histogram(dangle_1,32,'FaceColor','b','Normalization','probability')
histogram(dangle_2,32,'FaceColor','r','Normalization','probability')
xlim([-1,1]*5)
xlabel('flow angle - boundary angle (°)'); ylabel('probability')

% close all
% 
% figure
% hold on
% % plot(x2*scl.x,rad2deg(angle(x2)),'k')
% plot(bounds_x2*scl.x,dangle_A,'k')
% plot(bounds_x2*scl.x,dangle_1,'b')
% plot(bounds_x2*scl.x,dangle_2,'r')
% grid
% legend('Brute force','unmasked harmonic fit','masked harmonic fit')
% xlabel(lbl.x); ylabel('flow angle - boundary angle (°)')

%%
close all

figure
set(gcf,'PaperPosition',[0,0,6.5,2]*2)
tlo = tiledlayout(1,3);
xlim_zoom = [60,90];
ylim_zoom = [-40,-30];
vs = 2;

nexttile
hold on
quiver(X2*scl.x,X3*scl.x,v2_A*scl.v*vs,v3_A*scl.v*vs,0,'.-r')
plot(bounds_x2*scl.x,bounds_x3(1,:)*scl.x,'b')
xlim(xlim_zoom); ylim(ylim_zoom)
pbaspect(ar)

nexttile
hold on
quiver(X2*scl.x,X3*scl.x,v2_1*scl.v*vs,v3_1*scl.v*vs,0,'.-r')
plot(bounds_x2*scl.x,bounds_x3(1,:)*scl.x,'b')
xlim(xlim_zoom); ylim(ylim_zoom)
pbaspect(ar)

nexttile
hold on
quiver(X2*scl.x,X3*scl.x,v2_2*scl.v*vs,v3_2*scl.v*vs,0,'.-r')
plot(bounds_x2*scl.x,bounds_x3(1,:)*scl.x,'b')
xlim(xlim_zoom); ylim(ylim_zoom)
pbaspect(ar)

%% happy?
% happy = input('Happy? (y/n) ','s');
% if strcmp(happy,'y')
%     fprintf("I'm glad you're happy.\n")

%     direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\isinglass_05_A';
%     phi = phitopA;
%     save(fullfile(direc,'ext','potential_map.mat'),'phi','E2_bg','E3_bg','dec')
    
%     direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\isinglass_05_B';
%     phi = phitopB;
%     save(fullfile(direc,'ext','potential_map.mat'),'phi','E2_bg','E3_bg','ns') 
    
%     direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\isinglass_06_nosc_noro';
%     phi = phitopC;
%     save(fullfile(direc,'ext','potential_map.mat'),'phi','E2_bg','E3_bg')
% else
%     fprintf('aww\n')
% end

%% gather up isinglass flow data
% load('../../public_html/Gemini3D/isinglass_05/latlonggeoandgeomfootpoint.mat')
% load('../../public_html/Gemini3D/isinglass_05/fused_flow.mat')
% latf = griddedInterpolant(trajTime,Lat100km);
% lonf = griddedInterpolant(trajTime,Long100km);
% mlatf = griddedInterpolant(trajTime,mLat);
% mlonf = griddedInterpolant(trajTime,mLong);
% lat_100km = latf(time);
% lon_100km = lonf(time);
% mlat = mlatf(time);
% mlon = mlonf(time);
% datetimes = datetime(2017,3,2,7,50,0)+seconds(time);
% v_geog_east = ew;
% v_geog_north = ns;
%
% RE = 6370e3;
% dt = 1e-3; % seconds
%
% r0 = (RE + 100e3)*ones(size(lat_100km)); % 100 km footpoints, yeesh
% theta0 = deg2rad(90-lat_100km);
% phi0 = deg2rad(lon_100km);
%
% dr = 0;
% dtheta = -v_geog_north*dt./r0;
% dphi = v_geog_east*dt./(r0.*sin(theta0));
%
% r1 = r0 + dr;
% theta1 = theta0 + dtheta;
% phi1 = phi0 + dphi;
%
% [U0,E0,N0] = gemini3d.geog2UEN(r0,rad2deg(phi0),90-rad2deg(theta0), ...
%     mean(theta0),mean(phi0));
% [U1,E1,N1] = gemini3d.geog2UEN(r1,rad2deg(phi1),90-rad2deg(theta1), ...
%     mean(theta0),mean(phi0));
% v_geom_east = (E1-E0)/dt;
% v_geom_north = (N1-N0)/dt;
%
% figure
% quiver(E0,N0,v_geom_east,v_geom_north)
% pbaspect([1,1,1])
%
% figure
% hold on
% plot(v_geom_east.^2+v_geom_north.^2)
% plot(v_geog_east.^2+v_geog_north.^2)
%
% theta = deg2rad(28);
% v_geom_east_r = cos(theta)*ew-sin(theta)*ns;
% v_geom_north_r = sin(theta)*ew+cos(theta)*ns;
%
% comment = "JVI + KAL best understanding as of september 2023. " + ...
%     "Flow vectors rotated from geographic to geomagnetic according " + ...
%     "to gemini3d.geog2UEN using dt=1e-3 seconds. " + ...
%     "This amounts to a geographic up axis rotatation of 28.0 - 28.4 °.";
% note = "gemini3d.goeg2UEN is a simple west+south rotation, not " + ...
%     "proper apex transform. Hardwired to Alaska, 2017?";
%
% save('../../public_html/Gemini3D/isinglass_05/flow_data.mat', ...
%     "time", "datetimes", ...
%     "mlat", "mlon", ...
%     "lat_100km", "lon_100km", ...
%     "v_geog_east", "v_geog_north", ...
%     "v_geom_east", "v_geom_north", ...
%     "comment", "note" ...
%     )
