do_dat = 1;
run_direc = "\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D";
runs = [...
    "aurora_E0-dependance_01_A"...
    ,"aurora_E0-dependance_01_B"...
    ,"aurora_null_02"...
    ,"aurora_E0-dependance_01_C"...
    ,"aurora_E0-dependance_01_D"...
    ,"aurora_E0-dependance_01_E"...
    ,"aurora_highQlowE0_01"...
    ,"aurora_sharc_02"...
    ,"aurora_E2BG_05"
    ];
lruns = length(runs);
alt_ref = 300;
lon_ref = 800;

xg = gemini3d.read.grid(fullfile(run_direc,runs(1)));

%%
ltubes = 11;
p0 = zeros(ltubes,3);
p0(:,1) = lon_ref;
p0(:,2) = linspace(20,180,ltubes);
p0(:,3) = alt_ref;
r0 = ones(1,ltubes)*200;
r1 = ones(1,ltubes)*20;

%%
dats = struct;
E0 = zeros(lruns,1);
t = 35990;
for i = 1:lruns
    direc = fullfile(run_direc,runs(i));
    cfg = gemini3d.read.config(direc);
    time = datetime([cfg.ymd,0,0,t]);
    if do_dat
        dats.(char(64+i)) = gemini3d.read.frame(direc,'time',time);
    end
    E0(i) = cfg.E_amp_l;
end

%%
p1 = zeros(lruns,ltubes);
p3 = zeros(lruns,ltubes);
d1 = zeros(lruns,ltubes);
d2 = zeros(lruns,ltubes);
d3 = zeros(lruns,ltubes);
for i = 1:lruns
    disp(runs(i))
    dat = dats.(char(64+i));
    if i == 7
        offset = [-lon_ref,0,0];
    else
        offset = [0,0,0];
    end
    for n = 1:ltubes
        lastwarn('','')
        tube = fluxtube(xg,dat,alt_ref,p0(n,:)+offset,r0(n),r1(n),reverse=1,res=64);
        [msg,id] = lastwarn();
        if not(isempty(msg))
            p3(i,n) = missing;
            d1(i,n) = missing;
            d2(i,n) = missing;
            d3(i,n) = missing;
        else
            c0 = tube.caps.start;
            c1 = tube.caps.end;
            verts = tube.vertices;
            vertmins = cell2mat(cellfun(@(v) min(v(:,3)),verts,'UniformOutput',false));
            cls0 = min(vertmins);
            cls1 = max(vertmins);
            ctr0 = mean(c0);
            ctr1 = mean(c1);
            p1(i,n) = mean([cls0,cls1]);
            p3(i,n) = ctr0(2);
            d1(i,n) = cls1 - cls0;
            d2(i,n) = ctr1(1)-ctr0(1);
            d3(i,n) = ctr1(2)-ctr0(2);
        end
    end
end

%%
colors = [zeros(1,lruns-3);linspace(0,1,lruns-3);zeros(1,lruns-3)]';
fts = 17;
ftn = 'consolas';

reset(0)
set(0,'defaultFigurePaperUnits','inches')
set(0,'defaultTiledlayoutPadding','tight')
set(0,'defaultTiledlayoutTileSpacing','tight')
setall(0,'FontName',ftn)
setall(0,'FontSize',fts)
setall(0,'Multiplier',1)
set(0,'defaultAxesFontSizeMode','manual')
set(0,'defaultSurfaceEdgeColor','flat')

folder = 'morphology';
fn = 'E0-dependancy_morp';

legend_strings = [...
    num2str(E0(1:end-3)/1e3)...
    ,repmat(' keV  ',lruns-3,1)...
    ,num2str(E0(1:end-3)/1e3)...
    ,repmat(' mW/m^2',lruns-3,1)...
    ];

figure
set(gcf,'PaperPosition',[0,0,18,4.5])
tlo = tiledlayout(1,4);

nexttile
hold on
for i=1:lruns-3
    plot(-d2(i,:),-d3(i,:),'-o','Color',colors(i,:))
end
plot(-d2(end-2,:),-d3(end-2,:),'-ob')
plot(-d2(end-1,:),-d3(end-1,:),'-or')
plot(-d2(end-0,:),-d3(end-0,:),'-om')
xlabel('E-W Deflection, \DeltaX_E [km]')
ylabel('N-S Deflection, \DeltaX_N [km]')
grid on

nexttile
hold on
for i=1:lruns-3
    plot(p3(i,:),d2(i,:)./d3(i,:),'-o','Color',colors(i,:))
end
plot(p3(end-2,:),d2(end-2,:)./d3(end-2,:),'-ob')
plot(p3(end-1,:),d2(end-1,:)./d3(end-1,:),'-or')
plot(p3(end,:),d2(end,:)./d3(end,:),'-om')
xlabel('FAC Inflection Distance [km]')
ylabel('Deflection Ratio, \DeltaX_E/\DeltaX_N')
grid on

nexttile
hold on
a = p3;
b = p1;
for i=[1,lruns-3]%1:lruns-3
    plot(a(i,:),b(i,:),'-o','Color',colors(i,:))
end
plot(a(end-2,:),b(end-2,:),'-ob')
plot(a(end-1,:),b(end-1,:),'-or')
plot(a(end-0,:),b(end-0,:),'-om')
xlabel('FAC Inflection Distance [km]')
ylabel('Avg. Closure Altitude [km]')
grid on

nexttile
hold on
a = p3;
b = d1;
for i=1:lruns-3
    plot(a(i,:),b(i,:),'-o','Color',colors(i,:))
end
plot(a(end-2,:),b(end-2,:),'-ob')
plot(a(end-1,:),b(end-1,:),'-or')
plot(a(end-0,:),b(end-0,:),'-om')
xlabel('FAC Inflection Distance [km]')
ylabel('Closure Alt. Range, \DeltaX_U [km]')
grid on
legend([legend_strings;" 1 keV 100 mW/m^2";" Sharc";" E_E = 10 mV/m"],'Location','northeast','FontSize',0.8*fts,'box','off')

if ~exist(fullfile('plots',folder),'dir')
    mkdir(fullfile('plots',folder));
end
filename = fullfile('plots',folder,[fn,'.png']);
disp(['Saving: ',filename])
saveas(gcf,filename)
close all

%%
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
