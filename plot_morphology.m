do_dat = 1;
run_direc = "\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D";
runs = [...
    "aurora_E0-dependance_01_A"...
    ,"aurora_E0-dependance_01_B"...
    ,"aurora_null_02"...
    ,"aurora_E0-dependance_01_C"...
    ,"aurora_E0-dependance_01_D"...
    ,"aurora_E0-dependance_01_E"...
    ,"aurora_sharc_02"...
    ,"aurora_highQlowE0_01"
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
p3 = zeros(lruns,ltubes);
d2 = zeros(lruns,ltubes);
d3 = zeros(lruns,ltubes);
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
for i = 1:lruns
    disp(runs(i))
    dat = dats.(char(64+i));
    if i == lruns
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
            d2(i,n) = missing;
            d3(i,n) = missing;
        else
            c0 = tube.caps.start;
            c1 = tube.caps.end;
            ctr0 = mean(c0);
            ctr1 = mean(c1);
            p3(i,n) = ctr0(2);
            d2(i,n) = ctr1(1)-ctr0(1);
            d3(i,n) = ctr1(2)-ctr0(2);
        end
    end
end

%%
colors = [zeros(1,lruns-2);linspace(0,1,lruns-2);zeros(1,lruns-2)]';
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

% pp = ["\Deltax_N","\Deltax_E"];
% legend_strings = strings(2*lruns,1);
% for p = 1:length(pp)
%     for i = 1:lruns
%         legend_strings(2*i-1+p-1) = pp(p) + " " + num2str(E0(i)/1e3) + " keV";
%     end
% end
legend_strings = [...
    num2str(E0(1:end-2)/1e3)...
    ,repmat(' keV  ',lruns-2,1)...
    ,num2str(E0(1:end-2)/1e3)...
    ,repmat(' mW/m^2',lruns-2,1)...
    ];

figure
set(gcf,'PaperPosition',[0,0,15,5])
tlo = tiledlayout(1,3);

nexttile
hold on
for i=1:lruns-2
    plot(p3(i,:),d2(i,:),'-o','Color',colors(i,:))
end
plot(p3(end-1,:),d2(end-1,:),'-ob')
plot(p3(end,:),d2(end,:),'-or')
title('East-West Deflection')
xlabel('FAC Inflection Distance [km]')
ylabel('Relative Distance, \DeltaX_E [km]')
legend([legend_strings;" 3 keV   3 mw/m^2";" 1 keV 100 mW/m^2"],'Location','southwest','FontSize',0.6*fts)
grid on

nexttile
hold on
for i=1:lruns-2
    plot(p3(i,:),d3(i,:),'-o','Color',colors(i,:))
end
plot(p3(end-1,:),d3(end-1,:),'-ob')
plot(p3(end,:),d3(end,:),'-or')
title('North-South Deflection')
xlabel('FAC Inflection Distance [km]')
ylabel('Relative Distance, \DeltaX_N [km]')
legend([legend_strings;" 3 keV   3 mw/m^2";" 1 keV 100 mW/m^2"],'Location','southwest','FontSize',0.6*fts)
grid on

nexttile
hold on
for i=1:lruns-2
    plot(p3(i,:),d2(i,:)./d3(i,:),'-o','Color',colors(i,:))
end
plot(p3(end-1,:),d2(end-1,:)./d3(end-1,:),'-ob')
plot(p3(end,:),d2(end,:)./d3(end,:),'-or')
title('Deflection Ratio')
xlabel('FAC Inflection Distance [km]')
ylabel('Deflection Ratio, \DeltaX_N/\DeltaX_E')
legend([legend_strings;" 3 keV   3 mw/m^2";" 1 keV 100 mW/m^2"],'Location','northwest','FontSize',0.6*fts)
grid on

% nexttile
% hold on
% for i=1:lruns-1
%     plot(d2(i,:),d3(i,:),'-o','Color',colors(i,:))
% end
% plot(d2(end,:),d3(end,:),'-or')
% title('Deflection Ratio')
% xlabel('FAC Inflection Distance [km]')
% ylabel('Deflection Ratio, \DeltaX_N/\DeltaX_E')
% legend([legend_strings;" 1 keV, 100 mW/m^2"],'Location','southeast','FontSize',0.6*fts)
% grid on

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


% figure
% set(gcf,'PaperPosition',[0,0,5,3])
% title(plot_title,'FontSize',fts,'FontName',ftn,'FontWeight','bold','Interpreter','none')
% subtitle(['alt = ',num2str(alt_rac_p),' ',units.x,', mlon = ',num2str(mlon_rac_p),'°'])
%
% hold on
% yyaxis left
% plot(x3_p,E3_p)
% plot(x3_p,-j1_p(alt_rid,:))
% yyaxis right
% plot(x3_p,SIGP_p)
% plot(x3_p,SIGH_p)
% hold off
% xlabel(north_label)
% legend(...
%     ['E_N [',units.e,']']...
%     ,['j_{||} [',units.j,']']...
%     ,['\Sigma_P [',units.S,']']...
%     ,['\Sigma_H [',units.S,']']...
%     )
% xlim(x3_range_p)
% grid on
%
% if ~exist(fullfile(direc,'plots',folder),'dir')
%     mkdir(direc,fullfile('plots',folder));
% end
% filename = fullfile(direc,'plots',folder,[filename_prefix,'_',suffix,'.png']);
% disp(['Saving: ',filename])
% saveas(gcf,filename)