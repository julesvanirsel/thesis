load('data\replicate_data.mat')
v_bg = [0,-500];

%%
MLAT = 90-squeeze(xg.theta(end,:,:))*180/pi;
MLON = squeeze(xg.phi(end,:,:))*180/pi;
x2 = double(xg.x2(3:end-2))';
x3 = double(xg.x3(3:end-2))';
[X2,X3] = ndgrid(x2,x3);
mlon_to_x2 = griddedInterpolant(MLON(:,1),x2);
mlat_to_x3 = griddedInterpolant(MLAT(1,:),x3);
x2_to_mlon = griddedInterpolant(x2,MLON(:,1));
x3_to_mlat = griddedInterpolant(x3,MLAT(1,:));
Bmag = abs(mean(xg.Bmag,'all'));

in_situ_A = in_situ;
[~,v2_int_A,v3_int_A] = tools.replicate(in_situ_A,image,xg,flow_bg=v_bg);
fv2 = griddedInterpolant(X2,X3,v2_int_A);
fv3 = griddedInterpolant(X2,X3,v3_int_A);

%%
close all

x2_A = mlon_to_x2(in_situ_A.pos(:,1));
x3_A = mlat_to_x3(in_situ_A.pos(:,2));
v2_A = in_situ_A.flow(:,1);
v3_A = in_situ_A.flow(:,2);

dx2 = 80e3;
% dx2 = 20;
dx3 = -10e3;
psi = deg2rad(-15);
noise_amp = 50;
x2_tmp = x2_A + dx2;
x3_tmp = x3_A + dx3;
x2_B = cos(psi)*x2_tmp - sin(psi)*x3_tmp;
x3_B = sin(psi)*x2_tmp + cos(psi)*x3_tmp;
x2_B = -x2_B;
v2_B = fv2(x2_B,x3_B) + v_bg(1) + noise_amp*randn(size(x2_B));
v3_B = fv3(x2_B,x3_B) + v_bg(2) + noise_amp*randn(size(x2_B));

in_situ_B.pos(:,1) = x2_to_mlon(x2_B);
in_situ_B.pos(:,2) = x3_to_mlat(x3_B);
in_situ_B.flow(:,1) = v2_B;
in_situ_B.flow(:,2) = v3_B;

figure(1)
hold on
pcolor(mlon_to_x2(image.pos(:,:,1)),mlat_to_x3(image.pos(:,:,2)),image.flux)
quiver(x2_A,x3_A,v2_A,v3_A,0,'.-r')
quiver(x2_B,x3_B,v2_B,v3_B,0,'.-g')

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
traj_A = griddedInterpolant(tools.minsmooth(x3_A(sort_ids_A)),x2_A);
traj_B = griddedInterpolant(tools.minsmooth(x3_B(sort_ids_B)),x2_B);

% x_A = traj_A(x3);
% x_B = traj_B(x3);
% width = x_B - x_A;
% slope_rel = 0.5;
% slope = slope_rel*repmat(width,[length(x2),1])/2;
% weight_A = (1-tanh((x2'-x_A-width/2)./slope))/2;
% weight_B = 1-weight_A;

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

slope = 5e-5*0+1;
weight_A = (1 + tanh(slope*(min_dist_B-min_dist_A)))/2;
weight_B = 1 - weight_A;

% [~,sort_ids_A] = sort(x2_A);
% [~,sort_ids_B] = sort(x2_B);
% traj_A = griddedInterpolant(tools.minsmooth(x2_A(sort_ids_A)),x3_A);
% traj_B = griddedInterpolant(tools.minsmooth(x2_B(sort_ids_B)),x3_B);
% 
% y_A = traj_A(x2)';
% y_B = traj_B(x2)';
% width = y_B - y_A;
% slope_rel = 0.5;
% slope = slope_rel*repmat(width,[1,length(x3)])/2;
% weight_A = (1-tanh((x3-y_A-width/2)./slope))/2;
% weight_B = 1-weight_A;

figure
hold on
pcolor(X2,X3,weight_A)
plot(x2_A,x3_A)
plot(x2_B,x3_B)
shading flat
colorbar
clim([0,1])

%%
[~,v2_int_B,v3_int_B] = tools.replicate(in_situ_B,image,xg,flow_bg=v_bg);

%%
close all

v2_int = v2_int_A.*weight_A + v2_int_B.*weight_B;
v3_int = v3_int_A.*weight_A + v3_int_B.*weight_B;

lim.x = [min(x2),max(x2)];
lim.y = [min(x3),max(x3)];

figure
tiledlayout(3,2)

nexttile
hold on
pcolor(X2,X3,v2_int_A+v_bg(1))
quiver(x2_A,x3_A,v2_A-v_bg(1),v3_A-v_bg(2),1,'.-r')
shading flat
colorbar
xlim(lim.x); ylim(lim.y)

nexttile
hold on
pcolor(X2,X3,v3_int_A+v_bg(2))
quiver(x2_A,x3_A,v2_A-v_bg(1),v3_A-v_bg(2),1,'.-r')
shading flat
colorbar
xlim(lim.x); ylim(lim.y)

nexttile
hold on
pcolor(X2,X3,v2_int_B+v_bg(1))
quiver(x2_B,x3_B,v2_B-v_bg(1),v3_B-v_bg(2),1,'.-r')
shading flat
colorbar
xlim(lim.x); ylim(lim.y)

nexttile
hold on
pcolor(X2,X3,v3_int_B+v_bg(2))
quiver(x2_B,x3_B,v2_B-v_bg(1),v3_B-v_bg(2),1,'.-r')
shading flat
colorbar
xlim(lim.x); ylim(lim.y)

nexttile
hold on
pcolor(X2,X3,v2_int+v_bg(1))
quiver(x2_A,x3_A,v2_A-v_bg(1),v3_A-v_bg(2),1,'.-r')
quiver(x2_B,x3_B,v2_B-v_bg(1),v3_B-v_bg(2),1,'.-r')
shading flat
colorbar
xlim(lim.x); ylim(lim.y)

nexttile
hold on
pcolor(X2,X3,v3_int+v_bg(2))
quiver(x2_A,x3_A,v2_A-v_bg(1),v3_A-v_bg(2),1,'.-r')
quiver(x2_B,x3_B,v2_B-v_bg(1),v3_B-v_bg(2),1,'.-r')
shading flat
colorbar
xlim(lim.x); ylim(lim.y)