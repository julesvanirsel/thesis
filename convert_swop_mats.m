load('data\swop\skymap.mat','magnetic_footpointing')
asi_glon = magnetic_footpointing.('110km').lon;
asi_glat = magnetic_footpointing.('110km').lat;
xlims = [min(asi_glon(:)),max(asi_glon(:))];
ylims = [min(asi_glat(:)),max(asi_glat(:))];
clear('magnetic_footpointing')

mat_dir = fullfile('data','swop');
mat_filenames = {dir(fullfile(mat_dir, 'SW*.mat')).name};
colors = [["630","428","558"];["r","b","g"];["L13","L15","L14"]];
colorcet = @jules.tools.colorcet;
num = 1;

for f = mat_filenames
    fprintf('%i / %i\n',num,length(mat_filenames)); num = num + 1;
    mat_file = fullfile(mat_dir,cell2mat(f));
    skip = false;
    fid0 = fopen(fullfile('data','swop','keep.txt'),'r');
    fid1 = fopen(fullfile('data','swop','discard.txt'),'r');
    while ~feof(fid0)
        l = fgetl(fid0);
        if l==-1
            continue
        end
        if mat_file == l(1:48)
            fprintf(' %s...keep\n',mat_file)
            skip = true;
        end
    end
    while ~feof(fid1)
        l = fgetl(fid1);
        if mat_file == l
            fprintf(' %s...discard\n',mat_file)
            skip = true;
        end
    end
    fclose all;
    if skip
        continue
    end
    dt = datetime('today') - datetime(mat_file(30:37),'InputFormat','uuuuMMdd');

    % load swarm data
    load(mat_file,'UTtime','Footpoint_Long','Footpoint_Lat','Flow_y','FAC_micro')
    UTtime = pad(string(num2str(UTtime')),10,'left','0');
    track_times = datetime(UTtime,'InputFormat','HHmmss.SSS') - dt;
    track_glon = Footpoint_Long;
    track_glat = Footpoint_Lat;
    track_ve = Flow_y * sign(mean(diff(track_glat)));
    % track_be = By * sign(mean(diff(track_glat)));
    % track_be = track_be - mean(track_be);
    track_fac = FAC_micro;
    clear('Footpoint_Long','Footpoint_Lat','Flow_y','FAC_micro')

    % load imager data
    load(mat_file,'Imagery')
    for c = colors
        asi_times.(c(2)) = datetime(cell2mat(fieldnames(Imagery.(c(1)))) ...
            ,'InputFormat','HHmmss')' - dt;
        asi_frames.(c(2)) = reshape(struct2array(Imagery.(c(1))), ...
            [size(asi_glon),length(asi_times.(c(2)))]);
    end
    clear('Imagery')

    % find approx. center time
    [~,track_id] = min(abs(track_glat - mean(asi_glat(:),'omitnan')));
    sim_time = track_times(track_id);
    for c = colors
        [~,asi_id.(c(2))] = min(abs(sim_time - asi_times.(c(2))));
        sim_frames = asi_frames.(c(2));
        sim_frame.(c(2)) = double(sim_frames(:,:,asi_id.(c(2))))/1000;
    end

    % plot
    close all
    figure('Position',[0,50,1540,840])
    tiledlayout(2,3,'TileSpacing','tight')
    sgtitle(sprintf('FILE = %s   TIME = %s',mat_file,sim_time),'Interpreter','none','FontSize',15)

    for c = colors
        frame = sim_frame.(c(2));
        time = asi_times.(c(2));
        times = time(asi_id.(c(2))+(-5:5));
        times.Format = 'HH:mm:ss';
        [~,track_ids] = min(abs(track_times - times));
        clims = [quantile(frame(:),0.01),quantile(frame(:),0.99)];
        % clims = [2350,2550];

        nexttile
        hold on
        pcolor(asi_glon,asi_glat,frame)
        plot(track_glon,track_glat,'m','LineWidth',0.6)
        scatter(track_glon(track_ids),track_glat(track_ids),60,'mx','LineWidth',0.6)
        for i = 1:2:11
            id = track_ids(i);
            text(track_glon(id)+0.5,track_glat(id),num2str(i),'Color','m','FontSize',9)
        end
        clb = colorbar;
        colormap(gca,colorcet(char(c(3))))
        clb.Label.String = sprintf('%s nm @ %s',c(1),times(6));
        clim(clims)
        xlim(xlims)
        ylim(ylims)
        xlabel('Geo. lon.')
        ylabel('Geo. lat.')
        pbaspect([1,1,1])
        set(gca,'FontSize',15)
    end
    nexttile
    hold on
    pcolor(asi_glon,asi_glat,sim_frame.g)
    quiver(track_glon,track_glat,track_ve*2e-3,zeros(size(track_ve)) ...
        ,0,'.-m','LineWidth',0.6)
    quiver(-157.5,69.3,2,0,0,'.-m','LineWidth',0.6)
    text(-155.3,69.3,'1 km/s','FontSize',10)
    clb = colorbar;
    colormap(gca,colorcet(char(c(3))))
    clb.Label.String = sprintf('%s nm @ %s',c(1),times(6));
    clim(clims)
    xlim(xlims)
    ylim(ylims)
    xlabel('Geo. lon.')
    ylabel('Geo. lat.')
    pbaspect([1,1,1])
    set(gca,'FontSize',15)

    nexttile
    hold on
    pcolor(asi_glon,asi_glat,sim_frame.g)
    quiver(track_glon,track_glat,track_fac,zeros(size(track_ve)) ...
        ,0,'.-c','LineWidth',0.6)
    quiver(-157.5,69.3,2,0,0,'.-c','LineWidth',0.6)
    text(-155.3,69.3,'2 uA/m^2','FontSize',10)
    clb = colorbar;
    colormap(gca,colorcet(char(c(3))))
    clb.Label.String = sprintf('%s nm @ %s',c(1),times(6));
    clim(clims)
    xlim(xlims)
    ylim(ylims)
    xlabel('Geo. lon.')
    ylabel('Geo. lat.')
    pbaspect([1,1,1])
    set(gca,'FontSize',15)

    nexttile
    hold on
    pcolor(asi_glon,asi_glat,sim_frame.g)
    plot(track_glon,track_glat,'w','LineWidth',0.6)
    plot(track_glon + track_fac,track_glat,'c','LineWidth',0.6)
    plot(track_glon + track_ve*2e-3,track_glat,'m','LineWidth',0.6)
    clb = colorbar;
    colormap(gca,colorcet(char(c(3))))
    clb.Label.String = sprintf('%s nm @ %s',c(1),times(6));
    clim(clims)
    xlim(xlims)
    ylim(ylims)
    xlabel('Geo. lon.')
    ylabel('Geo. lat.')
    pbaspect([1,1,1])
    set(gca,'FontSize',15)

    keep = false;
    while true
        keep_q = input('Keep? (y/n): ','s');
        if keep_q == 'y'
            keep = true;
            break
        elseif keep_q == 'n'
            break
        end
    end
    if keep
        while true
            grade = input('Grade? (1-3): ');
            time_id = input('Frame number? (1-11): ');
            if any(grade==1:3) && any(time_id==1:11)
                break
            end
        end
        fid = fopen(fullfile('data','swop','keep.txt'),'a+');
        fprintf(fid,'%s\t%s\t%s\n',mat_file,track_times(track_ids(time_id)),num2str(grade));
        saveas(gcf,[mat_file(1:end-3),'png'])
    else
        fid = fopen(fullfile('data','swop','discard.txt'),'a+');
        fprintf(fid,'%s\n',mat_file);
    end
    fclose all;
end

%%
clear
% swop_dir = 'D:Files\research\thesis\data\swop';
% swop_dir = fullfile('data','swop');
fid = fopen(fullfile('data','swop','keep.txt'),'r');
colors = [["630","428","558"];["r","b","g"];["L13","L15","L14"]];
colorcet = @jules.tools.colorcet;

load('data\swop\skymap.mat','magnetic_footpointing')
asi_glon = magnetic_footpointing.('110km').lon;
asi_glat = magnetic_footpointing.('110km').lat;
xlims = [min(asi_glon(:)),max(asi_glon(:))];
ylims = [min(asi_glat(:)),max(asi_glat(:))];
clear('magnetic_footpointing')

ok_mats = repmat(' ',24,71);
i = 1;
while ~feof(fid)
    l = fgetl(fid);
    ok_mats(i,:) = sprintf('%s',l);
    i = i+1;
end
[good_dates,sort_ids] = sort(datetime(ok_mats(:,50:end-2)));
ok_mats = ok_mats(sort_ids,:);
good_ids = ok_mats(:,end) ~= '1';
ok_mats = string(ok_mats);
good_mats = ok_mats(good_ids)';

for f = good_mats(1)
    f_data = strsplit(f,'\t');
    mat_file = char(f_data(1));
    time = datetime(f_data(2));
    dur = 60;

    load(mat_file,'UTtime','Footpoint_Long','Footpoint_Lat','Flow_xh','Flow_xv','Flow_y','FAC_micro')
    UTtime = pad(string(num2str(UTtime')),10,'left','0');
    dt = datetime('today') - datetime(mat_file(30:37),'InputFormat','uuuuMMdd');
    track_times = datetime(UTtime,'InputFormat','HHmmss.SSS')' - dt;
    track_glon = Footpoint_Long;
    track_glat = Footpoint_Lat;
    track_ve = double(Flow_y * sign(mean(diff(track_glat))));
    track_vnh = double(Flow_xh * sign(mean(diff(track_glat))));
    track_vnv = double(Flow_xv * sign(mean(diff(track_glat))));
    track_fac = FAC_micro;
    [~,track_ids(2)] = min(abs(time-track_times));
    [~,track_ids(1)] = min(abs(time-seconds(dur/2)-track_times));
    [~,track_ids(3)] = min(abs(time+seconds(dur/2)-track_times));
    clear('UTtime','Footpoint_Long','Footpoint_Lat','Flow_xh','Flow_xv','Flow_y','FAC_micro')

    load(mat_file,'Imagery')
    for c = colors
        asi_times = datetime(cell2mat(fieldnames(Imagery.(c(1)))) ...
            ,'InputFormat','HHmmss')' - dt;
        asi_frames = reshape(struct2array(Imagery.(c(1))), ...
            [size(asi_glon),length(asi_times)]);
        [~,asi_id] = min(abs(time-asi_times));
        asi_frame.(c(2)) = asi_frames(:,:,asi_id);
        asi_time.(c(2)) = asi_times(asi_id);
    end
    clear('Imagery')

    img(:,:,1) = double(asi_frame.r)';
    img(:,:,2) = double(asi_frame.g)';
    img(:,:,3) = double(asi_frame.b)';
    lb = quantile(img(:),0.01);
    ub = quantile(img(:),0.99);
    img = flipud((img - lb)/(ub - lb));

    clims.g = [quantile(asi_frame.g(:),0.01),quantile(asi_frame.g(:),0.99)];
    tlims = track_times(track_ids([1,3]))+[-1,1]*seconds(1);
    vlims = [-1,1]*2;

    time.Format = 'hh:mm:ss   MMM dd';
    sat = mat_file(22);

    close all
    figure('Position',[0,50,1540,840])
    tiledlayout(5,5,'TileSpacing','tight')
    sgtitle(sprintf('FILE = %s    TIME = %s    SAT = %s',mat_file,time,sat) ...
        ,'Interpreter','none','FontSize',15)

    nexttile(1,[3,2])
    imshow(img)

    nexttile(3,[3,3])
    hold on
    pcolor(asi_glon,asi_glat,asi_frame.g)
    plot(track_glon,track_glat,'w','LineWidth',0.6)
    plot(track_glon + track_fac,track_glat,'c','LineWidth',0.6)
    plot(track_glon + track_ve*2e-3,track_glat,'m','LineWidth',0.6)
    scatter(track_glon(track_ids(1)),track_glat(track_ids(1)),100,'gx','LineWidth',0.6)
    scatter(track_glon(track_ids(2)),track_glat(track_ids(2)),100,'wx','LineWidth',0.6)
    scatter(track_glon(track_ids(3)),track_glat(track_ids(3)),100,'yx','LineWidth',0.6)
    clb = colorbar;
    colormap(gca,colorcet(char(c(3))))
    clb.Label.String = sprintf('%s nm @ %s',c(1),asi_time.g);
    clim(clims.g)
    xlim(xlims)
    ylim(sort(track_glat(track_ids([1,3])))+[-1,1]*0.5)
    xlabel('Geo. lon.')
    ylabel('Geo. lat.')
    set(gca,'FontSize',15)

    nexttile(16,[1,5])
    hold on
    plot(track_times,track_ve/1e3,'m','LineWidth',0.6)
    plot(track_times,track_vnh/1e3,'r','LineWidth',0.6)
    plot(track_times,track_vnv/1e3,'--r','LineWidth',0.6)
    plot([track_times(track_ids(1)),track_times(track_ids(1))],[-1,1]*100,'--g')
    plot([track_times(track_ids(2)),track_times(track_ids(2))],[-1,1]*100,'--w')
    plot([track_times(track_ids(3)),track_times(track_ids(3))],[-1,1]*100,'--y')
    xlim(tlims)
    ylim(vlims)
    grid on
    xlabel('Swarm UT')
    ylabel('Flow (km/s)')
    set(gca,'FontSize',15,'Color',[1,1,1]*0,'GridColor','w','GridAlpha',0.6)
    legend(["v_E","v_{NH}","v_{NV}"],'Color','w','Location','Eastoutside')


    nexttile(21,[1,5])
    hold on
    plot(track_times,track_fac,'c','LineWidth',0.6)
    plot(track_times,track_ve/1e3,'m','LineWidth',0.6)
    plot([track_times(track_ids(1)),track_times(track_ids(1))],[-1,1]*100,'--g')
    plot([track_times(track_ids(2)),track_times(track_ids(2))],[-1,1]*100,'--w')
    plot([track_times(track_ids(3)),track_times(track_ids(3))],[-1,1]*100,'--y')
    xlim(tlims)
    ylim(vlims)
    grid on
    xlabel('Swarm UT')
    ylabel('FAC (uA/m^2)')
    set(gca,'FontSize',15,'Color',[1,1,1]*0,'GridColor','w','GridAlpha',0.6)
    legend(["FAC","v_E"],'Color','w','Location','Eastoutside')

    
    set(gcf, 'InvertHardCopy', 'off'); 
    time.Format = 'MMdd';
    saveas(gcf,fullfile('plots','swop',sprintf('swop_%s_%s.png',time,sat)))
end

close all

