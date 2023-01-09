runs = struct2cell(dir('./runs/aurora_*'));
runs = runs(1,:);
lruns = length(runs);

time = datetime(2015,2,1,0,0,36000);
folder = 'D:\Files\research\thesis\runs\';
direc = fullfile(folder,runs{1});
xg = gemini3d.read.grid(direc);
alt = xg.x1(3:end-2)/1e3;

%%
doplots = 1;

clm = colorcet('D1A');
data = nan(lruns,4);
for run = 1:lruns
    disp(['run: ',num2str(run)])
    direc = fullfile(folder,runs{run});

    cfg = gemini3d.read.config(direc);
    fid = fopen(cfg.nml);
    while not(feof(fid))
        s = fgetl(fid);
        if contains(s,'PhiWBG')
            Q = str2double(s(10:end-17));
        elseif contains(s,'W0BG')
            E0 = str2double(s(8:end-17))/1e3;
        end
    end
    
    data(run,1:2) = [E0,Q];

    try
        dat = gemini3d.read.frame(direc,'time',time,'var','ne');
        ne = squeeze(dat.ne(:,round(end/2)));
        [~,i] = max(ne);
        data(run,3:4) = [alt(i),log10(ne(i))];
    catch
        disp(['SKIPPED: ',direc])
    end

    if doplots && run>20 && run<31
        cl = clm(round(((run-20)/10)*255+1),:);
        figure(1)
        hold on
        plot(log10(ne),alt,'color',cl)
        ylim([80,300])
        yline(alt(i),'color',cl)
    end
    if doplots && run>30
        cl = clm(round(((run-30)/10)*255+1),:);
        figure(2)
        hold on
        plot(log10(ne),alt,'color',cl)
        ylim([80,300])
        yline(alt(i),'color',cl)
    end
end

%% get Q from ne at max alt
Qf = fit(data(1:20,4),log10(data(1:20,2)),'poly1');
nemax = linspace(10,13,10);

hold on
scatter(data(1:10,4),log10(data(1:10,2)),'r')
scatter(data(11:20,4),log10(data(11:20,2)),'b')
plot(nemax,Qf(nemax),'k')
xlabel('log_{10} n_{e,max} [m^{-3}]')
ylabel('log_{10} Q [mW/m^2]')
xlim([10,13])
ylim([-2,3])
legend( ...
    ['E_0 = ', num2str(data(1,1)), ' keV'] ...
    ,['E_0 = ', num2str(data(11,1)), ' keV'] ...
    ,'Location','northwest'...
    )

%% get E0 from Q and max alt
hold on
scatter(log10(data(21:30,1)),data(21:30,3),'r')
scatter(log10(data(31:40,1)),data(31:40,3),'b')
legend( ...
    ['Q = ', num2str(data(21,2)), ' keV'] ...
    ,['Q = ', num2str(data(31,2)), ' keV'] ...
    ,'Location','northwest'...
    )

%%
% direc=['//Dartfs-hpc/rc/lab/L/LynchK/public_html/pubdata/2021_clayton_jgr/data/Simulations/c5' ...
% direc = './20170302_28365.000000.dat';                % folder where model output data is located
% direc = './runs/c5';
% time = datetime([2017,03,02,0,0,28365]);        % what time is of interest
% xg = gemini3d.read.grid(direc);   % from mat_gemini, reads in the grid info

%%
% simdat = gemini3d.read.frame(direc, "time", time); % reads in the model data
