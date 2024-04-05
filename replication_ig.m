%%
direc = '//dartfs-hpc/rc/lab/L/LynchK/public_html/Gemini3D/isinglass_78';
cfg = gemini3d.read.config(direc);
direc_rep = fullfile(direc,'/ext');
cfg_rep = gemini3d.read.config(direc_rep);
xg_rep = gemini3d.grid.cartesian(cfg_rep);
load('data\replicate_data_isinglass.mat','in_situ','image')

%%
[phi,~,E2_bg,E3_bg] = jules.tools.replicate(in_situ,image,xg_rep ...
    ,flow_smoothing_window = 16 ...
    ,boundary_smoothing_window = 64 ...
    ,show_plots = false ...
    ,save_plots = [0 0 0] ...
    ,save_data = false ...
    ,direc = 'plots\paper0' ...
    ,suffix = '__' ...
    ,starting_letter = 'A' ...
    ,add_phi_background = false ...
    ,fit_harmonic = true ...
    ,num_replications = 512 ... %32 or 512
    ,arc_definition = "conductance" ...
    ,edge_method = "contour" ...
    ,do_rotate = 0 ...
    ,do_scale = 0 ...
    ,harmonic_mask = [8,8,20]*1e3 ...
    ,contour_values = [10.5,19.1] ...
    );
phi_old = phi;

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
%% *** old ***

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

%
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
