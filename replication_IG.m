%% init
scl.x = 1e-3;
scl.v = 1e-3;
unt.x = 'km';
unt.v = 'km/s';
unt.q = 'mW/m^2';
lim.x = [-1,1]*120;
lim.y = [-1,1]*60;
lim.v = [-1,1]*1.3;
lim.verr = lim.v/2;
clm.v = colorcet('D2');
clm.divv = colorcet('CBTD1');
clm.Q = colorcet('L19');
lbl.x = ['Magnetic east (',unt.x,')'];
lbl.y = ['Magnetic north (',unt.x,')'];
lbl.vx = ['Eastward flow (',unt.v,')'];
lbl.vy = ['Northward flow (',unt.v,')'];
sc1 = 10;
sc2 = 1.5;

plotting = false;
%#ok<*UNRCH>

%% unpack grid and configuration data
direc = '//dartfs-hpc/rc/lab/L/LynchK/public_html/Gemini3D/isinglass_05';
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

if plotting
    close all
    figure(1)
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
    colormap(clm.Q)
    xlim(lim.x)
    ylim(lim.y)
    xlabel(lbl.x)
    ylabel(lbl.y)
end

clear('glon','glat','mlon','mlat','Qit','E0it','outputdate')

%% replicate rocket data
it_isin = 15;
v2_bg = 0;
v3_bg = -500; % chi-by-eye constant background flow removal
load(fullfile(direc,'ext','flow_data.mat'))
x2_isin = mlon_to_x2(mlon);
x3_isin = mlat_to_x3(mlat);
v2_isin = v_geom_east - v2_bg;
v3_isin = v_geom_north - v3_bg;
datetime_isin = datetimes(it_isin);

fprintf('Particle data datetime: %s\n', datetime_part)
fprintf('Rocket data datetime: %s\n', datetime_isin)

% 0 = original, 1 = replicated, a = primary bound, b = secondary bound
traj0 = griddedInterpolant(flip(x2_isin),flip(x3_isin));
x0a = fzero(@(x)(traj0(x)-bound_a(x)),0);
x0b = fzero(@(x)(traj0(x)-bound_b(x)),0);
y0a = traj0(x0a);
y0b = traj0(x0b);
width0 = sqrt((x0b-x0a)^2+(y0b-y0a)^2);
beta = atan2(x0b-x0a,y0b-y0a); % angle b/w bound-traj inters. and vertical

num_reps = 256;
dxs = linspace(-120e3,190e3,num_reps); % eastward displacements
x2_isin_rep = nan(length(dxs),length(x2_isin));
x3_isin_rep = nan(length(dxs),length(x3_isin));
v2_isin_rep = nan(length(dxs),length(v2_isin));
v3_isin_rep = nan(length(dxs),length(v3_isin));
for i = 1:length(dxs)
    % translate
    dx = dxs(i);
    dy = bound_a(x0a+dx)-y0a;

    x2_isin_tra = x2_isin + dx;
    x3_isin_tra = x3_isin + dy;

    % determine width at position 1
    traj1 = griddedInterpolant(flip(x2_isin_tra),flip(x3_isin_tra));
    x1a = fzero(@(x)(traj1(x)-bound_a(x)),0);
    x1b = fzero(@(x)(traj1(x)-bound_b(x)),0);
    y1a = traj1(x1a);
    y1b = traj1(x1b);
    width1 = sqrt((x1b-x1a)^2+(y1b-y1a)^2);

    % rotate about p1a by beta
    x2_isin_rot = cos(beta)*(x2_isin_tra - x1a) ...
        - sin(beta)*(x3_isin_tra - y1a) + x1a;
    x3_isin_rot = sin(beta)*(x2_isin_tra - x1a) ...
        + cos(beta)*(x3_isin_tra - y1a) + y1a;

    % scale about p1a
    scale = width1/width0;
    x2_isin_scl = x2_isin_rot;
    x3_isin_scl = scale*(x3_isin_rot - y1a) + y1a;

    % rotate back
    x2_isin_rep(i,:) = cos(-beta)*(x2_isin_scl - x1a) ...
        - sin(-beta)*(x3_isin_scl - y1a) + x1a;
    x3_isin_rep(i,:) = sin(-beta)*(x2_isin_scl - x1a) ...
        + cos(-beta)*(x3_isin_scl - y1a) + y1a;

    % rotate flows to be tangent to bound_a
    alpha = angle(x1a);
    v2_isin_rep(i,:) = cos(alpha)*v2_isin - sin(alpha)*v3_isin;
    v3_isin_rep(i,:) = sin(alpha)*v2_isin + cos(alpha)*v3_isin;

    if i==50 && plotting
        close all
        figure(1)
        hold on
        pcolor(X2_part*scl.x,X3_part*scl.x,Q_map_part); shading flat
        colormap("gray")
        plot(bounds_x2*scl.x,bounds_x3(1,:)*scl.x,'b')
        plot(bounds_x2*scl.x,bounds_x3(2,:)*scl.x,'b')
        scatter(x0a*scl.x,y0a*scl.x,50,'Filled','y')
        scatter(x0b*scl.x,y0b*scl.x,50,'Filled','y')
        scatter(x1a*scl.x,y1a*scl.x,50,'Filled','y')
        scatter(x1b*scl.x,y1b*scl.x,50,'Filled','y')
        quiver(x2_isin*scl.x,x3_isin*scl.x, ...
            v2_isin*scl.v*sc1,v3_isin*scl.v*sc1,0,'.-g')
        plot(x2_isin_tra*scl.x,x3_isin_tra*scl.x,'b')
        plot(x2_isin_rot*scl.x,x3_isin_rot*scl.x,'r')
        plot(x2_isin_scl*scl.x,x3_isin_scl*scl.x,'b--')
        quiver(x2_isin_rep(i,:)*scl.x,x3_isin_rep(i,:)*scl.x, ...
            v2_isin*scl.v*sc1,v3_isin*scl.v*sc1,0,'.-r')
        xlim(lim.x)
        ylim(lim.y)
        xlabel(lbl.x)
        ylabel(lbl.y)
        pbaspect([1,1,1])
    end
end

if plotting
    figure(2)
    hold on
    pcolor(X2_part*scl.x,X3_part*scl.x,Q_map_part); shading flat
    colormap("gray")
    plot(bounds_x2*scl.x,bounds_x3(1,:)*scl.x,'b')
    plot(bounds_x2*scl.x,bounds_x3(2,:)*scl.x,'b')
    quiver(x2_isin_rep*scl.x,x3_isin_rep*scl.x, ...
        v2_isin_rep*scl.v,v3_isin_rep*scl.v,0,'.-r')
    quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,0,'.-g')
    xlim(lim.x)
    ylim(lim.y)
    xlabel(lbl.x)
    ylabel(lbl.y)
    pbaspect([1,1,1])
end

clear('lon_100km','lat_100km','mlon','mlat', ...
    'v_geom_east','v_geom_north','v_geog_north','v_geog_north', ...
    'time','datetimes','comment','note')

%% interpolate replicated flow data
fv2 = scatteredInterpolant(x2_isin_rep(:),x3_isin_rep(:),v2_isin_rep(:));
fv3 = scatteredInterpolant(x2_isin_rep(:),x3_isin_rep(:),v3_isin_rep(:));
v2_map = fv2(X2,X3);
v3_map = fv3(X2,X3);
E2 =  v3_map*Bmag; % v = ExB/B^2
E3 = -v2_map*Bmag;
E2_bg =  v3_bg*Bmag;
E3_bg = -v2_bg*Bmag;

if plotting
    close all
    figure(1)
    hold on
    pcolor(X2*scl.x,X3*scl.x,v2_map*scl.v); shading flat
    colorbar; colormap(clm.v); clim(lim.v)
    quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v*sc1,v3_isin*scl.v*sc1,0,'.-k')
    xlim(lim.x)
    ylim(lim.y)
    xlabel(lbl.x)
    ylabel(lbl.y)
    pbaspect([1,1,1])

    figure(2)
    hold on
    pcolor(X2*scl.x,X3*scl.x,v3_map*scl.v); shading flat
    colorbar; colormap(clm.v); clim(lim.v)
    quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v*sc1,v3_isin*scl.v*sc1,0,'.-k')
    xlim(lim.x)
    ylim(lim.y)
    xlabel(lbl.x)
    ylabel(lbl.y)
    pbaspect([1,1,1])
end

%% generate potential map - option A
pseudo_basis = false;
usepar = false;
look_for_file = true;
suf = 'v7';
dec = [1,1];
numf = 64;

if pseudo_basis
    bag_x2 = X2(1:dec(1):end,1:dec(2):end);
    bag_x3 = X3(1:dec(1):end,1:dec(2):end);
    bag_v2 = v2_map(1:dec(1):end,1:dec(2):end);
    bag_v3 = v3_map(1:dec(1):end,1:dec(2):end);
    bag_x = [bag_x2(:),bag_x3(:)];
    bag_v = [bag_v2(:),bag_v3(:)];
    boundary = [bounds_x2;bounds_x3(1,:)]';

    if plotting
        close all

        figure(1)
        hold on
        pcolor(X2*scl.x,X3*scl.x,v2_map*scl.v); shading flat
        colorbar; colormap(clm.v); clim(lim.v)
        quiver(bag_x2*scl.x,bag_x3*scl.x,bag_v2*scl.v,bag_v3*scl.v,0,'.-k')

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

        phitopA = tools.basic_reconstruct(v2_map,v3_map,xg ...
            ,usepar=usepar,decimation=dec,ic_iteration=true);
        save(mat_fn,'phitopA','dec','usepar')
    end
end

%% generate potential map - option B
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
    phitopB = zeros([size(E2),prod(ns)]);
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
                    dot(E2(s2:d2:ix2,s3), dx2(s2:d2:ix2)*d2) ...
                    + dot(E3(ix2,s3:d3:ix3), dx3(s3:d3:ix3)*d3) ...
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

%% generate potential map - option C
load('forjewelz_threwee.mat')

%% plot potentials A and B
if plotting || true
    close all
    figure(1)
    tiledlayout(1,2)
    nexttile
    pcolor(X2,X3,phitopA); shading flat; colorbar
    nexttile
    pcolor(X2,X3,phitopB); shading flat; colorbar
end

%% plot interpolated flow against options A and B
add_bg = 1;
suff = 'v1';

if range(dx2) < 1e-3
    [E2_new_A,E3_new_A] = gradient(-phitopA',mean(dx2),mean(dx3));
    [E2_new_B,E3_new_B] = gradient(-phitopB',mean(dx2),mean(dx3));
else
    [E2_new_A,E3_new_A] = gradient(-phitopA',dx2,dx3);
    [E2_new_B,E3_new_B] = gradient(-phitopB',dx2,dx3);
end
v2_new_A = -E3_new_A'/Bmag;
v2_new_B = -E3_new_B'/Bmag;
v3_new_A =  E2_new_A'/Bmag;
v3_new_B =  E2_new_B'/Bmag;

dv2dx2 = diff(v2_map(:,2:end),1,1)./DX2(2:end,2:end);
dv3dx3 = diff(v3_map(2:end,:),1,2)./DX3(2:end,2:end);

close all
reset(0)
set(0,'defaultFigurePaperUnits','inches')
set(0,'defaultTiledlayoutPadding','tight')
set(0,'defaultTiledlayoutTileSpacing','tight')
set(0,'defaultSurfaceEdgeColor','flat')

figure(1)
set(gcf,'PaperPosition',[0,0,9,6.5])
tlo = tiledlayout(3,4);

% row 1
nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v2_map+add_bg*v2_bg)*scl.v)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
colormap(clm.v); clim(lim.v)
ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
xticks([])
title('Interpolated')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v2_new_A+add_bg*v2_bg)*scl.v)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
colormap(clm.v); clim(lim.v)
xlim(lim.x); ylim(lim.y)
xticks([]); yticks([])
title('Potential fit')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v2_new_B+add_bg*v2_bg)*scl.v)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
colormap(clm.v); clim(lim.v)
xlim(lim.x); ylim(lim.y)
xticks([]); yticks([])
title('Averaged path-integrated')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v2_map_sol+add_bg*v2_bg)*scl.v)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
colormap(clm.v); clim(lim.v)
clb = colorbar; clb.Label.String = lbl.vx;
xlim(lim.x); ylim(lim.y)
xticks([]); yticks([])
title('Helmholtz decomposition')

% row 2
nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_map+add_bg*v3_bg)*scl.v)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
colormap(clm.v); clim(lim.v)
ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
xticks([])

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_new_A+add_bg*v3_bg)*scl.v)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
colormap(clm.v); clim(lim.v)
xlim(lim.x); ylim(lim.y)
xticks([]); yticks([])

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_new_B+add_bg*v3_bg)*scl.v)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
colormap(clm.v); clim(lim.v)
xlim(lim.x); ylim(lim.y)
xlim(lim.x); ylim(lim.y)
xticks([]); yticks([])

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_map_sol+add_bg*v3_bg)*scl.v)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
colormap(clm.v); clim(lim.v)
clb = colorbar; clb.Label.String = lbl.vy;
xlim(lim.x); ylim(lim.y)
xticks([]); yticks([])

% row 3
nexttile
hold on
pcolor(X2(2:end,2:end)*scl.x,X3(2:end,2:end)*scl.x,dv2dx2+dv3dx3)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
clb = colorbar('southoutside'); clb.Label.String = '\nabla\cdot v (s^{-1})';
xlabel(lbl.x); ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
title('Divergence of Interp.')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_map-v3_new_A)*scl.v)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
colormap(clm.v); clim(lim.verr)
clb = colorbar('southoutside'); clb.Label.String = ['\Delta ',lbl.vy];
xlabel(lbl.x)
xlim(lim.x); ylim(lim.y)
yticks([])
title('Interp. − Pot. fit')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_new_A-v3_new_B)*scl.v)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
colormap(clm.v); clim(lim.verr)
clb = colorbar('southoutside'); clb.Label.String = ['\Delta ',lbl.vy];
xlabel(lbl.x)
xlim(lim.x); ylim(lim.y)
yticks([])
title('Pot. fit − Avg. path int.')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_new_A-v3_map_sol)*scl.v)
quiver(x2_isin*scl.x,x3_isin*scl.x,v2_isin*scl.v,v3_isin*scl.v,'.-r');
colormap(clm.v); clim(lim.verr)
clb = colorbar('southoutside'); clb.Label.String = ['\Delta ',lbl.vy];
xlabel(lbl.x)
xlim(lim.x); ylim(lim.y)
yticks([])
title('Pot. fit − Helmholtz dec.')

if pseudo_basis
    filename = sprintf('phitopA_%i_%i_%i_%s.png',dec(1),dec(2),numf,suf);
else
    filename = sprintf('phitopA_%i_%i_%s.png',dec(1),dec(2),suf);
end
filename = fullfile('plots',filename);
fprintf('Saving file:c %s\n',filename)
saveas(gcf,filename)

%% happy?
% load('data\phitopA_1_1_v7.mat')
% phi = phitopA;
% save(fullfile(direc,'ext','potential_map_01.mat'),'phi','E2_bg','E3_bg') 

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
% figure(1)
% quiver(E0,N0,v_geom_east,v_geom_north)
% pbaspect([1,1,1])
%
% figure(2)
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
