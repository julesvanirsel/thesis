%% generate simulation inputs
% nE0 = 10;
% nQ = 10;
% [E0g,Qg] = ndgrid(10.^linspace(2.001,4,nE0),10.^linspace(-1,2,nQ));
% 
% for i = 1:nE0
%     for j = 1:nQ
%         E0_tmp = E0g(i,j);
%         Q_tmp = Qg(i,j);
%         runname = ['pt_',num2str(E0_tmp,'%1.2e'),'_',num2str(Q_tmp,'%1.2e')];
%         run_series('runs\preciptest\aurora_preciptest_base\',runname,["W0BG","PhiWBG"],[E0_tmp;Q_tmp])
%     end
% end
% fclose all;

nE0 = 4;
nQ = 4;
[E0g,Qg] = ndgrid(10.^linspace(log10(500),log10(5e4),nE0),10.^linspace(-1,2,nQ));

for i = 1:nE0
    for j = 1:nQ
        E0_tmp = E0g(i,j);
        Q_tmp = Qg(i,j);
        runname = ['pt_',num2str(E0_tmp,'%1.2e'),'_',num2str(Q_tmp,'%1.2e'),'_2008'];
        run_series('runs\preciptest_fang2010\aurora_preciptest_base\',runname,["W0BG","PhiWBG","flag_fang"],[E0_tmp;Q_tmp;2008],run=true)
    end
end
fclose all;

for i = 1:nE0
    for j = 1:nQ
        E0_tmp = E0g(i,j);
        Q_tmp = Qg(i,j);
        runname = ['pt_',num2str(E0_tmp,'%1.2e'),'_',num2str(Q_tmp,'%1.2e'),'_2010'];
        run_series('runs\preciptest_fang2010\aurora_preciptest_base\',runname,["W0BG","PhiWBG","flag_fang"],[E0_tmp;Q_tmp;-1],run=true)
    end
end
fclose all;

%% make data matrix
runs = struct2cell(dir('.\runs\preciptest\pt_*'));
runs = runs(1,:);
lruns = length(runs);

time = datetime(2015,2,1,0,0,36000);
folder = 'D:\Files\research\thesis\runs\preciptest\';
direc = fullfile(folder,runs{1});
xg = gemini3d.read.grid(direc);
alt = xg.x1(3:end-2)/1e3;

clm = colorcet('D1A');
data = nan(lruns,4);
ne_data = nan(lruns,length(alt));
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
        ne_data(run,:) = log10(ne);
        [~,i] = max(ne);
        data(run,3:4) = [alt(i),log10(ne(i))];
    catch
        disp(['SKIPPED: ',direc])
    end
end

% order data set
i = 1;
data_ord = nan(size(data));
for E0 = unique(data(:,1))'
    ids = data(:,1)==E0;
    data_tmp = data(ids,:);
    [~,o] = sort(data_tmp(:,2));
    ld = length(data_tmp);
    data_ord(i:ld+i-1,:) = data_tmp(o,:);
    i = i+ld;
end
data = data_ord;

%%
ids = not(isnan(data(:,3))) & data(:,3)<250 & data(:,4)>10;

E0sf = fit([data(ids,3),data(ids,4)],log10(data(ids,1)),'poly33');
Qsf = fit([data(ids,3),data(ids,4)],log10(data(ids,2)),'poly22');

fs = 12;

figure(1)
plot(E0sf,[data(ids,3),data(ids,4)],log10(data(ids,1)))
alpha 0.5
set(gca,'FontSize',fs)
xlabel('Altitude of Peak Electron Density [km]')
ylabel('log_{10} Peak Electron Density [m^{-3}]')
zlabel('log_{10} Average Energy [keV]')

figure(2)
plot(Qsf,[data(ids,3),data(ids,4)],log10(data(ids,2)))
alpha 0.5
set(gca,'FontSize',fs)
xlabel('Altitude of Peak Electron Density [km]')
ylabel('log_{10} Peak Electron Density [m^{-3}]')
zlabel('log_{10} Total Energy Flux [mW/m^2]')

%%
clmi = round(linspace(1,256,10));
tiledlayout(2,5)

for i=1:10
    nexttile
    data_tmp = data(10*(i-1)+1:10*(i-1)+10,:);
    ne_data_tmp = ne_data(10*(i-1)+1:10*(i-1)+10,:);
    hold on
    for j=1:10
        plot(ne_data_tmp(j,:),alt,'Color','k')
    end
end
