% direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\aurora_sharc_02\';
direc = '\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\aurora_null_02\';
dat = gemini3d.read.frame(direc,'time',datetime(2015,2,1,0,0,36000));
xg = gemini3d.read.grid(direc);

%%
[X2,X3] = ndgrid(xg.x2(3:end-2)/1e3,xg.x3(3:end-2)/1e3);

%%
[~,precip_fn,~] = fileparts(dat.filename);
precip = dir(fullfile(direc,'*particles',[char(precip_fn),'.*']));
if isempty(precip)
    precip = dir(fullfile(direc,'*','*particles',[char(precip_fn),'.*']));
end
Q = h5read(fullfile(precip.folder,precip.name),'/Qp')/1e3; % Q defaults with mW/m^2 units
E0 = h5read(fullfile(precip.folder,precip.name),'/E0p'); % eV
jpar = -squeeze(dat.J1(end,:,:));

%%
clm.U = 'L19';
clm.c = 'L17';
clm.j = 'D1A';
units.U = 'mW/m^2';
units.c = 'keV';
units.j = 'uA/m^2';
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
set(gcf,'PaperPosition',[0,0,5,6])
tlo = tiledlayout(3,1);

nexttile
pcolor(X2,X3,Q*1e3)
title('Precip. Energy Flux')
xlabel('')
ylabel('North [km]')
colormap(gca,colorcet(clm.U))
clb = colorbar;
clb.Label.String = ['Q [',units.U,']'];
clb.Ruler.TickLabelFormat = clb_fmt;
clb.Ruler.Exponent = clb_exp;

nexttile
pcolor(X2,X3,E0/1e3)
title('Charactaristic Energy')
xlabel('')
ylabel('North [km]')
colormap(gca,colorcet(clm.c))
clb = colorbar;
clb.Label.String = ['E_0 [',units.c,']'];
clb.Ruler.TickLabelFormat = clb_fmt;
clb.Ruler.Exponent = clb_exp;

nexttile
pcolor(X2,X3,jpar*1e6)
title('Field Aligned Current')
xlabel('East [km]')
ylabel('North [km]')
colormap(gca,colorcet(clm.j))
clb = colorbar;
clb.Label.String = ['j_{||} [',units.j,']'];
clb.Ruler.TickLabelFormat = clb_fmt;
clb.Ruler.Exponent = clb_exp;
clim([-2,2])

saveas(gcf,'null_input_example.png')
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