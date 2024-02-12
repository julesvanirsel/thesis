if not(all(arrayfun(@exist,["in_situ","radar","image","xg","cfg"])))
    load('data\replicate_data_swop_03.mat')
end

%%
MLAT = 90-squeeze(xg.theta(end,:,:))*180/pi;
MLON = squeeze(xg.phi(end,:,:))*180/pi;
x2 = double(xg.x2(3:end-2))';
x3 = double(xg.x3(3:end-2))';
[X2,X3] = ndgrid(x2,x3);
lx2 = xg.lx(2); lx3 = xg.lx(3);
mlon_to_x2 = griddedInterpolant(MLON(:,1),x2);
mlat_to_x3 = griddedInterpolant(MLAT(1,:),x3);
x2_to_mlon = griddedInterpolant(x2,MLON(:,1));
x3_to_mlat = griddedInterpolant(x3,MLAT(1,:));
Bmag = abs(mean(xg.Bmag,'all'));

sig = 11;
in_situ_A = in_situ;
[~,~,~,E2_bg,E3_bg,v2_int_A,v3_int_A] = jules.tools.replicate(in_situ_A,image,xg ...
    ,flow_smoothing_window = 8 ...
    ,do_lowpass = true ...
    ,boundary_smoothing_window = 30 ...
    ,show_plots = true ...
    ,save_plots = false ...
    ,direc = 'plots\paper0' ...
    ,suffix = 'in_situ' ...
    ,add_phi_background = false ...
    ,fit_harmonic = true ...
    ,num_replications = 512 ...
    ,arc_definition = "conductance" ...
    ,edge_method = "contour" ...
    ,do_rotate = true ...
    ,do_scale = true ...
    ,contour_values = [1,1]*sig ...
    ,harmonic_mask = [1,0,1]*50e3 ...
    );

%%
v2_bg = -E3_bg/Bmag;
v3_bg =  E2_bg/Bmag;

in_situ_B = radar;
[~,~,~,~,~,v2_int_B,v3_int_B] = jules.tools.replicate(in_situ_B,image,xg ...
    ,flow_smoothing_window = 1 ...
    ,do_lowpass = true ...
    ,boundary_smoothing_window = 30 ...
    ,show_plots = true ...
    ,save_plots = false ...
    ,direc = 'plots\paper0' ...
    ,suffix = 'radar' ...
    ,add_phi_background = false ...
    ,fit_harmonic = true ...
    ,num_replications = 512 ...
    ,arc_definition = "conductance" ...
    ,edge_method = "contour" ...
    ,do_rotate = true ...
    ,do_scale = true ...
    ,contour_values = [1,1]*sig ...
    ,harmonic_mask = [1,0,1]*50e3 ...
    ,flow_bg = [v2_bg,v3_bg] ...
    );

% fv2 = griddedInterpolant(X2,X3,v2_int_A);
% fv3 = griddedInterpolant(X2,X3,v3_int_A);

%%
close all

x2_A = mlon_to_x2(in_situ_A.pos(:,1));
x3_A = mlat_to_x3(in_situ_A.pos(:,2));
v2_A = in_situ_A.flow(:,1);
v3_A = in_situ_A.flow(:,2);

x2_B = mlon_to_x2(in_situ_B.pos(:,1));
x3_B = mlat_to_x3(in_situ_B.pos(:,2));
v2_B = in_situ_B.flow(:,1);
v3_B = in_situ_B.flow(:,2);

% dx2 = 80e3;
% % dx2 = 20e3;
% dx3 = -10e3;
% psi = deg2rad(-15);
% noise_amp = 50;
% x2_tmp = x2_A + dx2;
% x3_tmp = x3_A + dx3;
% x2_B = cos(psi)*x2_tmp - sin(psi)*x3_tmp;
% x3_B = sin(psi)*x2_tmp + cos(psi)*x3_tmp;
% x2_B = -x2_B;
% v2_B = fv2(x2_B,x3_B) + v_bg(1) + noise_amp*randn(size(x2_B));
% v3_B = fv3(x2_B,x3_B) + v_bg(2) + noise_amp*randn(size(x2_B));
% 
% in_situ_B.pos(:,1) = x2_to_mlon(x2_B);
% in_situ_B.pos(:,2) = x3_to_mlat(x3_B);
% in_situ_B.flow(:,1) = v2_B;
% in_situ_B.flow(:,2) = v3_B;

sc = 1e2;
figure(1)
hold on
pcolor(mlon_to_x2(image.pos(:,:,1)),mlat_to_x3(image.pos(:,:,2)),image.flux)
quiver(x2_A,x3_A,v2_A*sc,v3_A*sc,0,'.-r')
quiver(x2_B,x3_B,v2_B*sc,v3_B*sc,0,'.-g')

figure(2)
hold on
plot(x3_A,v2_A,':r')
plot(x3_A,v3_A,':b')
plot(x3_B,v2_B,'r')
plot(x3_B,v3_B,'b')

%%
close all

x2_A = mlon_to_x2(in_situ_A.pos(:,1));
x3_A = mlat_to_x3(in_situ_A.pos(:,2));
x2_B = mlon_to_x2(in_situ_B.pos(:,1));
x3_B = mlat_to_x3(in_situ_B.pos(:,2));
lA = length(x2_A);
lB = length(x2_B);

[~,sort_ids_A] = sort(x3_A);
[~,sort_ids_B] = sort(x3_B);
traj_A = griddedInterpolant(jules.tools.minsmooth(x3_A(sort_ids_A)),x2_A);
traj_B = griddedInterpolant(jules.tools.minsmooth(x3_B(sort_ids_B)),x2_B);

distance_A = nan([size(X2),lA]);
for i=1:lA
    distance_A(:,:,i) = sqrt((X2-x2_A(i)).^2+(X3-x3_A(i)).^2);
end
distance_B = nan([size(X2),lB]);
for i=1:lB
    distance_B(:,:,i) = sqrt((X2-x2_B(i)).^2+(X3-x3_B(i)).^2);
end

min_dist_A = min(distance_A,[],3);
min_dist_B = min(distance_B,[],3);

slope = 1e-5;%*0+1;
weight_A = (1 + tanh(slope*(min_dist_B-min_dist_A)))/2;
weight_B = 1 - weight_A;
fprintf('max weight = %f\n',max(weight_A(:)))

[~,ids] = min(abs(weight_A-0.5));
x2_splt = smooth(x2(ids),1e3);
x3_splt = x3;

% weight_A = ones(size(weight_A));
% weight_B = zeros(size(weight_B));

figure
hold on
pcolor(X2,X3,weight_A)
plot(x2_A,x3_A)
plot(x2_B,x3_B)
plot(x2_splt,x3_splt,'k')
shading flat
colorbar
clim([0,1])

%%
% [~,v2_int_B,v3_int_B] = tools.replicate(in_situ_B,image,xg,flow_bg=v_bg);

%%
close all

v2_int = v2_int_A.*weight_A + v2_int_B.*weight_B;
v3_int = v3_int_A.*weight_A + v3_int_B.*weight_B;

scl.x = 1e-3; scl.v = 1e-3; scl.qv = 1e-1;
lim.x = [min(x2),max(x2)]*scl.x;
lim.y = [min(x3),max(x3)]*scl.x;
lim.v = [-1,1]*0.45;
clm.v = 'D2';

colorcet = @jules.tools.colorcet;

v_bg = [0,0];

v2_A_p = (v2_A-v2_bg)*scl.qv;
v3_A_p = (v3_A-v3_bg)*scl.qv;
v2_B_p = (v2_B-v2_bg)*scl.qv;
v3_B_p = (v3_B-v3_bg)*scl.qv;

divv = divergence(X3,X2,v3_int-v3_bg,v2_int-v2_bg);

close all
figure
set(gcf,'Position',[0,60,1920,1080-150]);
tiledlayout(4,2)

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,v2_int_A*scl.v)
quiver(x2_A*scl.x,x3_A*scl.x,v2_A_p,v3_A_p,0,'.-r')
shading flat
clb = colorbar;
clb.Label.String = 'v (km/s)';
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xlabel('M. east (km)'); ylabel('M. north (km)')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,v2_int_B*scl.v)
quiver(x2_B*scl.x,x3_B*scl.x,v2_B_p,v3_B_p,0,'.-r')
shading flat
clb = colorbar;
clb.Label.String = 'v (km/s)';
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xlabel('M. east (km)'); ylabel('M. north (km)')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,v3_int_A*scl.v)
quiver(x2_A*scl.x,x3_A*scl.x,v2_A_p,v3_A_p,0,'.-r')
shading flat
clb = colorbar;
clb.Label.String = 'v (km/s)';
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xlabel('M. east (km)'); ylabel('M. north (km)')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,v3_int_B*scl.v)
quiver(x2_B*scl.x,x3_B*scl.x,v2_B_p,v3_B_p,0,'.-r')
shading flat
clb = colorbar;
clb.Label.String = 'v (km/s)';
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xlabel('M. east (km)'); ylabel('M. north (km)')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,v2_int*scl.v)
quiver(x2_A*scl.x,x3_A*scl.x,v2_A_p,v3_A_p,0,'.-r')
quiver(x2_B*scl.x,x3_B*scl.x,v2_B_p,v3_B_p,0,'.-r')
plot(x2_splt*scl.x,x3_splt*scl.x,'k','LineWidth',0.5)
shading flat
clb = colorbar;
clb.Label.String = 'v (km/s)';
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xlabel('M. east (km)'); ylabel('M. north (km)')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,v3_int*scl.v)
quiver(x2_A*scl.x,x3_A*scl.x,v2_A_p,v3_A_p,0,'.-r')
quiver(x2_B*scl.x,x3_B*scl.x,v2_B_p,v3_B_p,0,'.-r')
plot(x2_splt*scl.x,x3_splt*scl.x,'k','LineWidth',0.5)
shading flat
clb = colorbar;
clb.Label.String = 'v (km/s)';
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xlabel('M. east (km)'); ylabel('M. north (km)')

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,divv*1e3)
quiver(x2_A*scl.x,x3_A*scl.x,v2_A_p,v3_A_p,0,'.-r')
quiver(x2_B*scl.x,x3_B*scl.x,v2_B_p,v3_B_p,0,'.-r')
plot(x2_splt*scl.x,x3_splt*scl.x,'k','LineWidth',0.5)
shading flat
clb = colorbar;
clb.Label.String = 'dv (mHz)';
% colormap(gca,colorcet(clm.dv))
xlim(lim.x); ylim(lim.y)
xlabel('M. east (km)'); ylabel('M. north (km)')

%% generate potential map
% translated from Alex Mule's python code Oct 9, 2023
close all
E2_int =  v3_int*Bmag; % v = ExB/B^2
E3_int = -v2_int*Bmag;

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

%% plot final flow fields
scl.x = 1e-3; scl.v = 1e-3; scl.dv = 1e3; scl.Q = 1e0; scl.p = 1e-3;
unt.x = 'km'; unt.v = 'km/s'; unt.dv = 'mHz'; unt.Q = 'mW/m^2'; unt.p = 'kV';
clm.v = 'D2'; clm.dv = 'CBD1'; clm.Q = 'L19'; clm.p = 'D10';
lim.x = [-1,1]*120; lim.y = [-1,1]*59;  lim.v = [-1,1]*1.3; lim.dv = [-1,1]*0.1;

lbl.x = sprintf('Magnetic east (%s)',unt.x);
lbl.y = sprintf('Magnetic north (%s)',unt.x);
lbl.vx = sprintf('Eastward flow (%s)',unt.v);
lbl.vy = sprintf('Northward flow (%s)',unt.v);

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

figure
set(gcf,'PaperPosition',[0,0,9,6.5])
tiledlayout(3,3);

qnt = 0.95;
max_v = quantile(abs([v2(:)+v_bg(1);v3(:)+v_bg(2)]),qnt);
max_dv = quantile(abs([divv(:);divv_int(:)]),qnt);
max_p = quantile(abs(phi(:)),qnt);
lim.v = [-1,1]*max_v*scl.v;
lim.dv = [-1,1]*max_dv*scl.dv;
lim.p = [-1,1]*max_p*scl.p;

% row 1
nexttile
title('Interpolated')
hold on
pcolor(X2*scl.x,X3*scl.x,(v2_int+v_bg(1))*scl.v)
quiver(x2_A*scl.x,x3_A*scl.x,v2_A*scl.v,v3_A*scl.v,'.-r')
quiver(x2_B*scl.x,x3_B*scl.x,v2_B*scl.v,v3_B*scl.v,'.-r')
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
quiver(x2_A*scl.x,x3_A*scl.x,v2_A*scl.v,v3_A*scl.v,'.-r')
quiver(x2_B*scl.x,x3_B*scl.x,v2_B*scl.v,v3_B*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xticks([]); yticks([])

nexttile
title('Difference & Potential map')
hold on
pcolor(X2*scl.x,X3*scl.x,v2_err*scl.v)
quiver(x2_A*scl.x,x3_A*scl.x,v2_A*scl.v,v3_A*scl.v,'.-r')
quiver(x2_B*scl.x,x3_B*scl.x,v2_B*scl.v,v3_B*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
clb = colorbar;
clb.Label.String = sprintf('\\Delta Eastward flow (%s)',unt.v);
xlim(lim.x); ylim(lim.y); clim(lim.v/3)
xticks([]); yticks([])

%row 2
nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,(v3_int+v_bg(2))*scl.v)
quiver(x2_A*scl.x,x3_A*scl.x,v2_A*scl.v,v3_A*scl.v,'.-r')
quiver(x2_B*scl.x,x3_B*scl.x,v2_B*scl.v,v3_B*scl.v,'.-r')
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
quiver(x2_A*scl.x,x3_A*scl.x,v2_A*scl.v,v3_A*scl.v,'.-r')
quiver(x2_B*scl.x,x3_B*scl.x,v2_B*scl.v,v3_B*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
xlim(lim.x); ylim(lim.y); clim(lim.v)
xticks([]); yticks([])

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,v3_err*scl.v)
quiver(x2_A*scl.x,x3_A*scl.x,v2_A*scl.v,v3_A*scl.v,'.-r')
quiver(x2_B*scl.x,x3_B*scl.x,v2_B*scl.v,v3_B*scl.v,'.-r')
colormap(gca,colorcet(clm.v))
clb = colorbar;
clb.Label.String = sprintf('\\Delta Northward flow (%s)',unt.v);
xlim(lim.x); ylim(lim.y); clim(lim.v/3)
xticks([]); yticks([])

% row 3
nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,divv_int*scl.dv)
colormap(gca,colorcet(clm.dv))
clb = colorbar;
clb.Label.String = sprintf('Divergence (%s)',unt.dv);
clb.Location = 'westoutside';
xlim(lim.x); ylim(lim.y); clim(lim.dv)
xlabel(lbl.x); ylabel(lbl.y)

nexttile
hold on
pcolor(X2*scl.x,X3*scl.x,divv*scl.dv)
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

