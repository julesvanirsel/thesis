run = 'runs/maeve_surge_02';
run = 'runs/maeve_drift_01';

if ~exist('xg','var')
    xg = gemini3d.read.grid(run);
end
x3 = xg.x3(3:end-2)/1e3;
lx1 = xg.lx(1);
lx2 = xg.lx(2);

cfg = gemini3d.read.config(run);
ymd = cfg.ymd;
UTsec0 = cfg.UTsec0;

t=150;
time = datetime([ymd,[0 0 UTsec0+t]]);
dat = gemini3d.read.frame(run,'Time',time);
[dir,fn,ext] = fileparts(dat.filename);
Q = h5read(dir+filesep+'OSSE_particles'+filesep+fn+ext,'/Qp');
E0 = h5read(dir+filesep+'OSSE_particles'+filesep+fn+ext,'/E0p');

p1 = lx1;
p2 = round(lx2/2);

J1 = squeeze(dat.J1(p1,p2,:))/1e-6;
v2 = squeeze(dat.v2(p1,p2,:))/1e3;
Q = squeeze(Q(p2,:))';
E0 = squeeze(E0(p2,:))'/1e3;

J1_max = round(max(J1));
v2_max = round(max(v2));
Q_max = round(max(Q));
E0_max = round(max(E0));

J1 = J1/J1_max;
v2 = v2/v2_max;
Q = Q/Q_max;
E0 = E0/E0_max;

%%
fz = 11;

set(gcf,'units','inches','OuterPosition',[1 1 5 4])

hold on
plot(x3,J1,'r')
plot(x3,Q,'b')
plot(x3,E0,'b--')
plot(x3,v2,'k')
hold off
xlim(400*[-1 1])
ylim(1.1*[min(min([J1 Q E0 v2])) max(max([J1 Q E0 v2]))])
xlabel('North [km]','FontSize',fz)
ylabel('Normalized Values','FontSize',fz)
set(gca,'FontSize',fz)
pbaspect([2 1 1])
grid on
legend(...
    'j_{||} / '+string(J1_max)+' uA/m^2'...
    ,'Q / '+string(Q_max)+' mW/m^2'...
    ,'E_0 / '+string(E0_max)+' keV'...
    ,'v_{east} / '+string(v2_max)+' km/s'...
    ,'fontsize',0.7*fz,'location','southeast'...
)

print(dir+filesep+fn+'_'+string(p1)+'_'+string(p2)+'_cut.png','-dpng');

