direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\aurora_sharc_02\';
dat = gemini3d.read.frame(direc,'time',datetime(2015,2,1,0,0,36000));
xg = gemini3d.read.grid(direc);

%%
[X1,X2,X3] = ndgrid(xg.x1(3:end-2)/1e3,xg.x2(3:end-2)/1e3,xg.x3(3:end-2)/1e3);

jpar = -squeeze(dat.J1(end,:,:));
v2 = dat.v2;
v3 = dat.v3;

buf = 1.05;
qnt=0.99;
j1_range = buf*[-1,1]*quantile(abs(jpar(:)),qnt);

[hsv_map_clb,hsv_mlon,hsv_mlon_map,hsv_alt,hsv_alt_map] =...
    tools.hsv_params(v2,v3,X3,X2,X1,300,0,700);

%%
clm.j = 'D1A';
units.j = 'uA/m^2';
units.v = 'km/s';
clb_fmt = '%+ 6.1f';
clb_exp = 0;
ftn = 'Consolas';
fts = 17;

reset(0)
set(0,'defaultFigurePaperUnits','inches')
set(0,'defaultTiledlayoutPadding','tight')
set(0,'defaultTiledlayoutTileSpacing','tight')
setall(0,'FontName',ftn)
setall(0,'FontSize',fts)
setall(0,'Multiplier',1)
set(0,'defaultAxesFontSizeMode','manual')
set(0,'defaultSurfaceEdgeColor','flat')

figure
set(gcf,'PaperPosition',[0,0,4,6])

pcolor(squeeze(X2(end,:,:)),squeeze(X3(end,:,:)),hsv_alt);
xlabel('East [km]')
ylabel('North [km]')
colormap(gca,hsv_alt_map)
clb = colorbar;
colormap(clb,hsv_map_clb)
clb.Limits = [0,1];
clb.Ticks = [0,1/4,1/2,3/4,1];
clb.TickLabels = {'W','S','E','N','W'};
clb.Label.String = ['Sat. at 0.7 ',units.v];
clim([0,1])
xlim([200-400,200+400])
ylim([-190,360])

saveas(gcf,'flowclosure.png')
close all

function setall(obj,property_suffix,val)
property_suffix = char(property_suffix);
l = length(property_suffix);
properties = fieldnames(get(obj,'factory'));
for n = 1:numel(properties)
    p = properties{n};
    if strcmp(p(end-l+1:end),property_suffix)
        set(obj,['default',p(8:end)],val)
    end
end
end