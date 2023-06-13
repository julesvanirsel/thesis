% rnm = 'null_02';
% rnm = 'sharc_02';
% rnm = 'bent_03';
rnm = 'E3BG_02';
% rnm = 'longsharc_01';

clear('dat')
clear('xg')

%%
direc = ['\\Dartfs-hpc\rc\lab\L\LynchK\public_html\Gemini3D\aurora_',rnm,'\'];
alt_ref = 300;
lon_ref = 600;
fs = 20;
xlims = [-1,1]*880;
ylims = [-190,360];
zlims = [80,alt_ref*1.05];
units.x = 'km';
units.j = 'uA/m^2'; clm.j = 'D1A';
units.n = 'm^{-3}'; clm.n = 'L9';
units.v = 'km/s';   clm.v = 'D2';
j1_range = [-1,1]*1.6;
ne_range = [9.9,11.5];
v2_range = [-1,1];
vv = [45,25];
% vv = [90,0];
iso_val = log10(2e11);
iso_alp = 0.3;

cfg = gemini3d.read.config(direc);
ymd = cfg.ymd;
UTsec0 = cfg.UTsec0;
tdur = cfg.tdur;
dtout = cfg.dtout;

if not(exist('xg','var'))
    xg = gemini3d.read.grid(direc);
end
MLAT = 90-squeeze(xg.theta(1,:,:))*180/pi;
MLON = squeeze(xg.phi(1,:,:))*180/pi;
x = double(xg.x2(3:end-2)/1e3);
y = double(xg.x3(3:end-2)/1e3);
z = double(xg.x1(3:end-2)/1e3);
[Xm,Ym,Zm] = meshgrid(x,y,z);
lz = length(z);

% load('../thesis/runs/c5_new/inputs/particles/particles.mat')
% [MLATp,MLONp] = ndgrid(mlat,mlon);

% for UTsec = UTsec0+5:dtout:UTsec0+tdur-15
for UTsec = UTsec0+150
    time = datetime(ymd) + seconds(UTsec);
%     time.Format = 'yyyyMMdd''T''HHmmss.SSS';
%     filename_prefix = [char(time),'UT'];
    if not(exist('dat','var'))
        dat = gemini3d.read.frame(direc,'time',time);
    end

    jz = permute(dat.J1,[2,3,1])*1e6;
    ne = permute(log10(dat.ne),[2,3,1]);
    v1 = permute(dat.v2*1e-3,[2,3,1]);
    v2 = permute(dat.v2*1e-3,[2,3,1]);
    v3 = permute(dat.v3*1e-3,[2,3,1]);
%     Qgi = griddedInterpolant(MLATp,MLONp,Qit(:,:,UTsec-UTsec0+1));
%     E0gi = griddedInterpolant(MLATp,MLONp,E0it(:,:,UTsec-UTsec0+1));
%     Q = Qgi(MLAT,MLON);
%     E0 = E0gi(MLAT,MLON)/1e3;

    v1 = permute(v1,[2,1,3]);
    v2 = permute(v2,[2,1,3]);
    v3 = permute(v3,[2,1,3]);
    
    if strcmp(rnm,'null_02')
        p0 = [[lon_ref,60,alt_ref];[lon_ref,150,alt_ref]];
        r0 = [1,1]*200;
        r1 = [1,1]*20;
        rev = [1,1];
        res = [1,1]*64;
        colors = [[1.0,0.7,0.0];[0.0,0.7,0.0]];
    elseif strcmp(rnm,'sharc_02')
        p0 = [[lon_ref,60,alt_ref];[lon_ref,160,alt_ref]];
        r0 = [1,1]*200;
        r1 = [1,1]*20;
        rev = [1,1];
        res = [1,1]*64;
        colors = [[1.0,0.7,0.0];[0.0,0.7,0.0]];
    elseif strcmp(rnm,'bent_03')
        p0 = [[lon_ref,60+50,alt_ref];[lon_ref,150+50,alt_ref]];
        r0 = [1,1]*200;
        r1 = [1,1]*20;
        rev = [1,1];
        res = [1,1]*64;
        colors = [[1.0,0.7,0.0];[0.0,0.7,0.0]];
    elseif strcmp(rnm,'E3BG_02')
        p0 = [[lon_ref,60,alt_ref];[lon_ref,150,alt_ref]];
        r0 = [1,1]*200;
        r1 = [1,1]*20;
        rev = [1,1];
        res = [1,1]*64;
        colors = [[1.0,0.7,0.0];[0.0,0.7,0.0]];
    end
    ntubes = size(p0,1);
    for n = 1:ntubes
        tubes.(char(64+n)) = fluxtube(xg,dat,alt_ref,p0(n,:),r0(n),r1(n)...
            ,reverse=rev(n),color=colors(n,:),res=res(n));
    end

    figure
    set(gcf,'PaperUnits','inches','PaperPosition',[0,0,11,6])
    %         title([runname,' at ',title_time,' UT'],'FontSize',fs*2,'FontWeight','bold','Interpreter','none')
    t = tiledlayout(1,1,'TileSpacing','compact');
    axj = axes(t); %#ok<LAXES>
    axn = axes(t); %#ok<LAXES>
    axt = axes(t); %#ok<LAXES>
    axa = [axj,axn,axt];

    set(axa,'FontSize',fs)
    set(axj,'XColor','none','YColor','none','ZColor','none','ZTick',alt_ref)
    set(axn,'Color','none','XColor','none','YColor','none','ZColor','none','XTick',lon_ref)
    set(axt,'Color','none','XGrid','on','YGrid','on','ZGrid','on')

    xlim(axj,xlims)
    xlim(axn,xlims + (lon_ref-xlims(1)))
    xlim(axt,xlims)
    ylim(axa,ylims)
    zlim(axj,zlims + (alt_ref-zlims(1)))
    zlim(axn,zlims)
    zlim(axt,zlims)
    view(axa,vv)
%     ar = [range(xlims),range(ylims),range(zlims)];
    ar = [250,128,88];
    pbaspect(axj,ar)
    pbaspect(axn,ar)
    pbaspect(axt,ar)
    set(axa,'YTick',[0,200])
    hold(axa,'on')

%     ida = 1+9*8;
%     vx = zeros(size(v1));
%     vy = zeros(size(v1));
%     vz = zeros(size(v1));
%     vx(:,:,ida) = v2(:,:,ida)*100;
%     vy(:,:,ida) = v3(:,:,ida)*100;
%     quiver3(axt,ds(Xm),ds(Ym),ds(Zm),ds(vx),ds(vy),ds(vz),0,'Color',[1,0,1],'LineWidth',1)
%     annotation('textarrow',[0.65+0.018,0.65],[0.85,0.85] ...
%         ,'String',' = 1 km/s','HeadLength',4,'HeadWidth',4,'LineWidth',1,'Color',[1,0,1],'FontSize',fs)
    
    slice(axj,Xm,Ym,Zm,permute(jz,[2,1,3]),[],[],alt_ref);
    colormap(axj,colorcet(clm.j))
    shading(axj,'flat')
    clim(axj,j1_range)
    clb = colorbar(axj);
    clb.Label.String = ['j_{||} [',units.j,']'];
    clb.FontSize = fs;
    clb.Position = [0.87,0.15,0.017,0.36];

%     slice(axj,Xm,Ym,Zm,repmat(Q,[1,1,lz]),[],[],alt_ref);
%     colormap(axj,colorcet('L1'))
%     shading(axj,'flat')
%     clim(axj,[0,42])
%     clb = colorbar(axj);
%     clb.Label.String = 'Q [mW/m^2]';
%     clb.FontSize = fs;
%     clb.Position = [0.86,0.15,0.017,0.36];

    slice(axn,Xm,Ym,Zm,permute(ne,[2,1,3]),lon_ref,[],[]);
    colormap(axn,colorcet(clm.n))
    shading(axn,'flat')
    clim(axn,ne_range)
    clb = colorbar(axn);
    clb.Label.String = ['log_{10} n_e [',units.n,']'];
    clb.FontSize = fs;
    clb.Position = [0.87,0.55,0.017,0.36];
    clb.Ticks = [10,11];

%     clri = round(255*(iso_val-ne_range(1))/range(ne_range))+1;
%     clm_tmp = colorcet(clm.n);
%     ne_iso = permute(ne,[2,1,3]);
%     s = isosurface(Xm,Ym,Zm,ne_iso,iso_val);
%     p = patch(s);
%     isonormals(Xm,Ym,Zm,ne_iso,p)
%     set(p,'FaceColor',clm_tmp(clri,:));
% %     set(p,'FaceColor',[0.6,0.3,0.8])
% %     set(p,'FaceColor',[1,0.5,0])
%     set(p,'FaceAlpha',iso_alp)
%     set(p,'EdgeColor','none');
% 
%     ctrs = contourslice(axn,Xm,Ym,Zm,permute(ne,[2,1,3]),lon_ref,[],[],[iso_val,13]);
%     set(ctrs,'EdgeColor','k','LineWidth',1)

    for n = 1:ntubes
        color = colors(n,:);
        tube = tubes.(char(64+n));
        verts = tube.vertices;
        c0 = tube.caps.start;
        c1 = tube.caps.end;
        in0 = tube.flux.area.in;
        in1 = tube.flux.area.out;

        shadow = nan(size(in0));
        shadow(in0) = 0;
        shadow(in1) = 0;
        shadow = repmat(shadow,[1,1,lz]);

        pl0 = plot3(axt,c0(:,1),c0(:,2),c0(:,3));
        pl0shd = plot3(axt,c0(:,1),c0(:,2),ones(size(c0))*zlims(1),'--');
        pl1 = plot3(axt,c1(:,1),c1(:,2),c1(:,3));
        pl1shd = plot3(axt,c1(:,1),c1(:,2),ones(size(c1))*zlims(1),'--');
        stl = streamline(axt,verts);
%         quiver3(c1(1:4:end,1),c1(1:4:end,2),c1(1:4:end,3)-1,zeros(res(n)/4,1),zeros(res(n)/4,1),1*ones(res(n)/4,1),0)
%         if any(not(isnan(shadow)),'all')
%             shd = slice(axj,Xm,Ym,Zm,permute(shadow,[2,1,3]),[],[],alt_ref);
%         end
%         plot3(axt,[1,1,1]*lon_ref,[ylims,ylims(2)],[zlims(1),zlims],'--','Color',[0,1,0]);
%         plot3(axt,[xlims(1),xlims],[ylims,ylims(2)],[1,1,1]*alt_ref,'--','Color',[0,0,1]);

        shading(axj,'flat')
        set(pl0,'Color',color,'LineWidth',2)
        set(pl0shd,'Color',color,'LineWidth',2)
        set(pl1,'Color',color,'LineWidth',2)
        set(pl1shd,'Color',color,'LineWidth',2)
%         set(shd,'FaceAlpha',0.5)
        set(stl,'Color',[color,0.5],'LineWidth',1)
    end

    xlabel(['East [',units.x,']'],'FontSize',fs)
    ylabel(['North [',units.x,']'],'FontSize',fs)
    zlabel(['Up [',units.x,']'],'FontSize',fs)

    fn = ['D:\Files\research\FINESST\',rnm,'.png'];
%     fn = 'D:\Files\research\FINESST\null.png';
    saveas(gcf,fn)
    disp(fn)
    close all
end

% function dr = ds(d)
% qs = 8;
% dr = d(1:qs:end,1:qs:end,1:qs:end);
% end