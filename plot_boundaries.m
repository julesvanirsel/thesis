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

load('data\flux_grad.mat')
bound2 = bound;
load('data\cond_cont.mat')
bound1 = bound;

figure
set(gcf,'PaperPosition',[0,0,13.2,4.8])
tiledlayout(1,2)
ltr = 65;

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8); ltr = ltr+1;
hold on
pcolor(X2_imag*scl.x,X3_imag*scl.x,arc.^ap*scl.arc)
plot(bound_pts*scl.x,bound1.A(bound_pts)*scl.x,'k')
plot(bound_pts*scl.x,bound1.B(bound_pts)*scl.x,'--k')
plot(bound_pts*scl.x,bound2.A(bound_pts)*scl.x,'Color',[0,0.9,0],'LineWidth',lw*0.8)
plot(bound_pts*scl.x,bound2.B(bound_pts)*scl.x,'--','Color',[0,0.9,0],'LineWidth',lw*0.8)
colormap(colorcet(clm.arc))
clb = colorbar;
clb.Label.String = 'Pedersen conductance (S)';
clb.Location = 'northoutside';
xlabel(lbl.x); ylabel(lbl.y)
xlim(lim.x); ylim(lim.y)
pbaspect(ar)

nexttile
text(0.04,0.9,char(ltr),'units','normalized','FontSize',fts*0.8);
hold on
pcolor(X2_imag(2:end-1,2:end-1)*scl.x,X3_imag(2:end-1,2:end-1)*scl.x,edges)
plot(bound_pts*scl.x,bound.A(bound_pts)*scl.x,'k')
plot(bound_pts*scl.x,bound.B(bound_pts)*scl.x,'--k')
plot(bound_pts*scl.x,bound2.A(bound_pts)*scl.x,'Color',[0,0.9,0],'LineWidth',lw)
plot(bound_pts*scl.x,bound2.B(bound_pts)*scl.x,'--','Color',[0,0.9,0],'LineWidth',lw)
colormap(colorcet(clm.arc))
clb = colorbar;
clb.Label.String = 'Sobel edges (a.u.)';
clb.Location = 'northoutside';
xlabel(lbl.x)
yticks([])
xlim(lim.x); ylim(lim.y)
pbaspect(ar)

filename = 'boundaries.png';
saveas(gcf,fullfile('plots','paper0',filename))
