%% load simulation
direc.A = "../../public_html/Gemini3D/isinglass_74/";
direc.B = "../../public_html/Gemini3D/isinglass_74_noscale_norotate/";
cfg.A = gemini3d.read.config(direc.A);
cfg.B = gemini3d.read.config(direc.B);
time = cfg.A.times(end);
xg = gemini3d.read.grid(direc.A);
dat.A = gemini3d.read.frame(direc.A,'time',time);
dat.B = gemini3d.read.frame(direc.B,'time',time);

%% unpack grid
MLAT = 90-squeeze(xg.theta)*180/pi;
MLON = squeeze(xg.phi)*180/pi;
x2 = xg.x2(3:end-2);
x3 = xg.x3(3:end-2);
[X2,X3] = ndgrid(x2,x3);
lx1 = xg.lx(1); lx2 = xg.lx(2); lx3 = xg.lx(3);

j1 = zeros([2,lx1,lx2,lx3]);
jA = zeros([2,lx2,lx3]);
jB = zeros([2,lx2,lx3]);
jC = zeros([2,lx2,lx3]);
for i = 1:2
    s = char(64+i);

    % formatting simulation data
    j1(i,:,:,:) = dat.(s).J1;

    phi = dat.(s).Phitop;
    [~,E2,E3] = gemscr.postprocess.pot2field(xg,phi);
    [~,~,SIGP,SIGH] = jules.tools.load_conductances(direc.(s),time,dat.(s),cfg.(s),xg);
    SIGH = -SIGH;
    
    % add background electric fields
    % E2 = E2 + cfg.(s).ap_ExBg;
    % E3 = E3 + cfg.(s).ap_EyBg;

    % E0_UTsecs = UTsec0 + (0:dtE0:tdur);
    % [~,E0_i] = min(abs(E0_UTsecs-UTsec));
    % E0_time = datetime(ymd) + seconds(E0_UTsecs(E0_i));
    E0_fn = fullfile(direc.A,cfg.A.E0_dir,[char(gemini3d.datelab(time)),'.h5']);
    E2_BG = mean(h5read(E0_fn,'/Exit'),'all');
    E3_BG = mean(h5read(E0_fn,'/Eyit'),'all');
    E2 = E2 + E2_BG;
    E3 = E3 + E3_BG;
    
    % implicit simulation data
    [~,dX2] = gradient(X2);
    [dX3,~] = gradient(X3);
    [~,dE22] = gradient(squeeze(E2(1,:,:)));
    [dE33,~] = gradient(squeeze(E3(1,:,:)));
    [dSIGP3,dSIGP2] = gradient(SIGP);
    [dSIGH3,dSIGH2] = gradient(SIGH);
    jA(i,:,:) = SIGP.*(dE22./dX2+dE33./dX3); % SIGP*div_perp(E)
    jB(i,:,:) = (dSIGP2.*squeeze(E2(1,:,:))./dX2 + dSIGP3.*squeeze(E3(1,:,:))./dX3); % grad(SIGP).E_perp
    jC(i,:,:) = (dSIGH2.*squeeze(E3(1,:,:))./dX2 - dSIGH3.*squeeze(E2(1,:,:))./dX3); % grad(SIGH).bxE_perp
end

%% plot
close all
fts = 20; % fontsize
ftn = 'Arial';
lw = 1.4;

reset(0)
set(0,'defaultFigurePaperUnits','inches')
set(0,'defaultTiledlayoutPadding','tight')
set(0,'defaultTiledlayoutTileSpacing','tight')
set(0,'defaultLineLineWidth',lw)
jules.tools.setall(0,'FontName',ftn)
jules.tools.setall(0,'FontSize',fts)
jules.tools.setall(0,'Multiplier',1)
set(0,'defaultAxesFontSizeMode','manual')
set(0,'defaultSurfaceEdgeColor','flat')

colorcet = @jules.tools.colorcet;

scl.x = 1e-3; scl.j = 1e6;
lim.x = [-1,1]*90; lim.y = [-1,1]*58;
units.j = 'uA/m^2';
clm.j = 'D1A';

j1_p = j1*scl.j;
jA_p = jA*scl.j; jB_p = jB*scl.j; jC_p = jC*scl.j;

lbl.x = 'Mag. E (km)';
lbl.y = 'Mag. N (km)';
j1_range_p = [-1,1]*90;

ar = [range(x2),range(x3),range(x3)];

load('data\boundaries.mat')
bound_x2 = linspace(lim.x(1),lim.x(2),256)*1e3;

figure(1)
set(gcf,'PaperPosition',[0,0,13.2,6.8])
tlo = tiledlayout(2,3);

% row 1
nexttile
title('div(E) term')
hold on
pcolor(X2/1e3,X3/1e3,squeeze(jA_p(2,:,:)))
plot(bound_x2/1e3,bound.A(bound_x2)/1e3-3,'--k')
plot(bound_x2/1e3,bound.B(bound_x2)/1e3,'--k')
ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
xticks([])
colormap(gca,colorcet(clm.j))
clim(j1_range_p*1.4)
pbaspect(ar)

nexttile
title('grad(\Sigma_P) term')
hold on
pcolor(X2/1e3,X3/1e3,squeeze(jB_p(2,:,:)))
plot(bound_x2/1e3,bound.A(bound_x2)/1e3-3,'--k')
plot(bound_x2/1e3,bound.B(bound_x2)/1e3,'--k')
xlim(lim.x); ylim(lim.y)
xticks([]); yticks([])
colormap(gca,colorcet(clm.j))
clim(j1_range_p*0.3)
pbaspect(ar)

nexttile
title('grad(\Sigma_H) term')
hold on
pcolor(X2/1e3,X3/1e3,squeeze(0*jA_p(2,:,:)+0*jB_p(2,:,:)+jC_p(2,:,:)))
plot(bound_x2/1e3,bound.A(bound_x2)/1e3-3,'--k')
plot(bound_x2/1e3,bound.B(bound_x2)/1e3,'--k')
xlim(lim.x); ylim(lim.y)
xticks([]); yticks([])
colormap(gca,colorcet(clm.j))
clim(j1_range_p*0.7)
pbaspect(ar)

% row 2
nexttile
hold on
pcolor(X2/1e3,X3/1e3,squeeze(jA_p(1,:,:)))
plot(bound_x2/1e3,bound.A(bound_x2)/1e3-3,'--k')
plot(bound_x2/1e3,bound.B(bound_x2)/1e3,'--k')
xlabel(lbl.x); ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
colormap(gca,colorcet(clm.j))
clim(j1_range_p*1.4)
clb = colorbar;
clb.Label.String = 'j_{||} (uA/m^2)';
clb.Location = 'southoutside';
pbaspect(ar)

nexttile
hold on
pcolor(X2/1e3,X3/1e3,squeeze(jB_p(1,:,:)))
plot(bound_x2/1e3,bound.A(bound_x2)/1e3-3,'--k')
plot(bound_x2/1e3,bound.B(bound_x2)/1e3,'--k')
xlabel(lbl.x)
xlim(lim.x); ylim(lim.y)
yticks([])
colormap(gca,colorcet(clm.j))
clim(j1_range_p*0.3)
clb = colorbar;
clb.Label.String = 'j_{||} (uA/m^2)';
clb.Location = 'southoutside';
pbaspect(ar)

nexttile
hold on
pcolor(X2/1e3,X3/1e3,squeeze(0*jA_p(1,:,:)+0*jB_p(1,:,:)+jC_p(1,:,:)))
plot(bound_x2/1e3,bound.A(bound_x2)/1e3-3,'--k')
plot(bound_x2/1e3,bound.B(bound_x2)/1e3,'--k')
xlabel(lbl.x)
xlim(lim.x); ylim(lim.y)
yticks([])
colormap(gca,colorcet(clm.j))
clim(j1_range_p*0.7)
clb = colorbar;
clb.Label.String = 'j_{||} (uA/m^2)';
clb.Location = 'southoutside';
pbaspect(ar)

% nexttile
% pcolor(X2/1e3,X3/1e3,squeeze(jA_p(2,:,:)-jA_p(1,:,:)))
% xlabel(lbl.x); ylabel(lbl.y)
% xlim(lim.x); ylim(lim.y)
% colormap(gca,colorcet(clm.j))
% clim(j1_range_p)
% pbaspect(ar)
% 
% nexttile
% pcolor(X2/1e3,X3/1e3,squeeze(jB_p(2,:,:)-jB_p(1,:,:)))
% xlabel(lbl.x)
% xlim(lim.x); ylim(lim.y)
% yticks([])
% colormap(gca,colorcet(clm.j))
% clim(j1_range_p)
% pbaspect(ar)
% 
% nexttile
% pcolor(X2/1e3,X3/1e3,squeeze(jC_p(2,:,:)-jC_p(1,:,:)))
% xlabel(lbl.x)
% xlim(lim.x); ylim(lim.y)
% yticks([])
% colormap(gca,colorcet(clm.j))
% clim(j1_range_p)
% clb = colorbar;
% clb.Label.String = '\Delta j_{||} (uA/m^2)';
% pbaspect(ar)

saveas(gcf,fullfile('plots','paper0','continuity-comparison.png'))
close all
