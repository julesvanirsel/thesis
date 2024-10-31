direc = fullfile('..','..','public_html','Gemini3D','swop_0314_AC_02/');
cfg = gemini3d.read.config(direc);
xg = gemini3d.read.grid(direc);
xg = jules.tools.shrink(xg);
time = cfg.times(end);
dat = gemini3d.read.frame(direc,"time",time,'vars',["J1","J2","J3"]);

sats = strsplit(direc,'_');
sats = sats{end-1};

%%
x1 = double(xg.x1(3:end-2));
x2 = double(xg.x2(3:end-2));
x3 = double(xg.x3(3:end-2));
[X1,X2,X3] = ndgrid(x1,x2,x3);
fFAC = griddedInterpolant(X1,X2,X3,-dat.J1*1e6);
mlon_to_x2 = griddedInterpolant(xg.mlon,x2);
mlat_to_x3 = griddedInterpolant(xg.mlat,x3);

%%
swrm.direc = fullfile('data','paper2','swarm_data');
swrm.fns = {dir(fullfile(swrm.direc,'SW_*FAC*.h5')).name};
for fn = swrm.fns
    sat = fn{1}(12);
    t0 = datetime(fn{1}(20:34),'InputFormat','uuuuMMdd''T''HHmmss');
    t1 = datetime(fn{1}(36:50),'InputFormat','uuuuMMdd''T''HHmmss');
    if (t0 <= time) && (time <= t1)
        tmp = fullfile(swrm.direc,fn{1});
        swrm.(sat).h5fn = tmp;
        swrm.(sat).time = datetime(h5read(tmp,'/Timestamp'),'ConvertFrom','posixtime');
        swrm.(sat).x1 = h5read(tmp,'/GeodeticAltitude');
        swrm.(sat).x2 = mlon_to_x2(h5read(tmp,'/MagneticLongitude'))'; % NEEDS TO BE GEODETIC MOVING FORWARD
        swrm.(sat).x3 = mlat_to_x3(h5read(tmp,'/MagneticLatitude'))';
        swrm.(sat).FAC = h5read(tmp,'/FAC');
        swrm.(sat).ids = (time-seconds(30) < swrm.(sat).time) & (swrm.(sat).time < time+seconds(30));
    end
end

%%
clr.A = [1,0,0]; clr.B = [0,1/2,0]; clr.C = [0,0,1];
offset = 10e3;
factor = 0.1;
close all
figure
hold on
for s = sats
    ids = swrm.(s).ids;
    gem.(s).FAC = fFAC(swrm.(s).x1(ids), swrm.(s).x2(ids)-200e3, swrm.(s).x3(ids));
    plot(swrm.(s).x3(ids), swrm.A.FAC(ids), ...
        'DisplayName', ['Swarm ',s], 'Color', clr.(s))
    plot(swrm.(s).x3(ids) + offset, gem.(s).FAC*factor, ...
        '--', 'DisplayName', ['Gemini ',s], 'Color', clr.(s))
end
xlim([min(x3),max(x3)])
legend
