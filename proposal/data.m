h5info('Proposal\SW_EXPT_EFIC_TIICT_20140422_0844T0849_fsmi.h5','/event').Datasets.Name;
img = h5read('Proposal\SW_EXPT_EFIC_TIICT_20140422_0844T0849_fsmi.h5','/event/Images');
img = img(:,:,33);
lat = h5read('Proposal\SW_EXPT_EFIC_TIICT_20140422_0844T0849_fsmi.h5','/event/Skymap/FULL_MAP_LATITUDE');
lat = lat(1:end-1,1:end-1,2);
lon = h5read('Proposal\SW_EXPT_EFIC_TIICT_20140422_0844T0849_fsmi.h5','/event/Skymap/FULL_MAP_LONGITUDE');
lon = lon(1:end-1,1:end-1,2);

satlat = h5read('Proposal\SW_EXPT_EFIC_TIICT_20140422_0844T0849_fsmi.h5','/event/Footpoint_Lat');
satlon = h5read('Proposal\SW_EXPT_EFIC_TIICT_20140422_0844T0849_fsmi.h5','/event/Footpoint_Long');

fac = h5read('Proposal\SW_EXPT_EFIC_TIICT_20140422_0844T0849_fsmi.h5','/event/FAC_micro');
flo = h5read('Proposal\SW_EXPT_EFIC_TIICT_20140422_0844T0849_fsmi.h5','/event/Unofficial_crosstrack_flow');
bri = h5read('Proposal\SW_EXPT_EFIC_TIICT_20140422_0844T0849_fsmi.h5','/event/Brightness_Vector');
zer = zeros(length(fac),1);

fz=19;

set(gcf,'units','inches','OuterPosition',[1 1 14 7])
xlims = [58 63];
xtic = [59.4 60 60.66];
pre = 184:209;
ret = 209:228;

sp1 = subplot(3,2,1);
hold on
plot(satlat,fac,'k')
plot(satlat(pre),fac(pre),'b','LineWidth',2)
plot(satlat(ret),fac(ret),'r','LineWidth',2)
hold off
grid on
xlim(xlims)
ylim([-10 10])
yticks([-10 -5 0 5 10])
ylabel('FAC [uA/m^2]')
set(gca,'fontsize',fz,'Xticklabel',[])

sp2 = subplot(3,2,3);
hold on
plot(satlat,flo/1e3,'k')
plot(satlat(pre),flo(pre)/1e3,'b','LineWidth',2)
plot(satlat(ret),flo(ret)/1e3,'r','LineWidth',2)
hold off
grid on
xlim(xlims)
ylim([-0.7 0.7])
yticks([-0.5 0 0.5])
ylabel('v_E [km/s]')
set(gca,'fontsize',fz,'Xticklabel',[])

sp3 = subplot(3,2,5);
hold on
plot(satlat(6:595),bri(:,1)/1e4,'k')
plot(satlat(pre),bri(pre-5,1)/1e4,'b','LineWidth',2)
plot(satlat(ret),bri(ret-5,1)/1e4,'r','LineWidth',2)
hold off
grid on
xlim(xlims)
ylim([0.3 1])
yticks([0.4 0.6 0.8 1])
xlabel('Latitude [°]')
ylabel('Brightness [ul]')
set(gca,'fontsize',fz)

sp4 = subplot(3,2,[2,4,6]);
hold on
pcolor(lon,lat,img);
quiver(satlon,satlat,satlon+1e10*fac,satlon,2,'-b','ShowArrowHead','off');
quiver(satlon,satlat+0.02,satlon+1e10*flo,satlon+0.02,1,'-r','ShowArrowHead','off');
hold off
shading flat
caxis([3000 5000])
colormap('Gray')
xlim([241 255])
ylim([55 65])
yticks(linspace(50,70,11))
set(gca,'fontsize',fz)
xlabel('Longitude [°]')
ylabel('Latitude [°]')
pbaspect([1 1 1])

% print('Proposal/data.png','-dpng')
