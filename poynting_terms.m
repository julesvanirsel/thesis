run = 'runs\maeve_4\';

if ~exist('xg','var')
    xg = gemini3d.read.grid(run);
end
x1 = xg.x1(4:end-3);
x2 = xg.x2(4:end-3);
x3 = xg.x3(4:end-3);
[X1,X2,X3] = ndgrid(x1,x2,x3);

cfg = gemini3d.read.config(run);
ymd = cfg.ymd;
UTsec0 = cfg.UTsec0;
tdur = cfg.tdur;
dt = cfg.dtout;

for t=150
    time = datetime([ymd,[0 0 UTsec0+t]]);
    dat = gemini3d.read.frame(run,"Time",time);
end

%%
phi = dat.Phitop;
phi = permute(repmat(phi,1,1,xg.lx(1)),[3 1 2]);

E1 = -der(phi,1,xg);
E2 = -der(phi,2,xg);
E3 = -der(phi,3,xg);

Bmag = unique(xg.Bmag);


v2 = dat.v2(2:end-1,2:end-1,2:end-1);
v3 = dat.v3(2:end-1,2:end-1,2:end-1);
v2_new =  E3./xg.Bmag(2:end-1,2:end-1,2:end-1);
v3_new = -E2./xg.Bmag(2:end-1,2:end-1,2:end-1);

%%

figure(1)
pcolor(squeeze(X2(end,:,:)),squeeze(X3(end,:,:)),squeeze(v2_new(end,:,:)))
colorbar
shading flat

figure(2)
pcolor(squeeze(X2(end,:,:)),squeeze(X3(end,:,:)),squeeze(v2_new(end,:,:)-v2(end,:,:)))
colorbar
shading flat

% figure(3)
% pcolor(squeeze(X2(end,:,:)),squeeze(X3(end,:,:)),squeeze(v3(end,:,:)))
% colorbar
% shading flat

% figure(4)
% quiver(squeeze(X2(end,:,:)),squeeze(X3(end,:,:)),squeeze(E2(end,:,:)),squeeze(E3(end,:,:)),3)

function v = der(u,d,xg)
    dim = [0 0 0];
    dim(d) = 1;
    
    x = xg.('x'+char(d));
    x = x(3:end-2);
    dx = x(2:end)-x(1:end-1);
    if d==1
        DX = permute(repmat(dx,1,xg.lx(2)-1,xg.lx(3)-1),[1 2 3]);
    elseif d==2
        DX = permute(repmat(dx,1,xg.lx(1)-1,xg.lx(3)-1),[2 1 3]);
    elseif d==3
        DX = permute(repmat(dx,1,xg.lx(1)-1,xg.lx(2)-1),[2 3 1]);
    end
    
    h1 = DX(1:end-1,1:end-1,1:end-1);
    h2 = DX(1+dim(1):end-1+dim(1)...
           ,1+dim(2):end-1+dim(2)...
           ,1+dim(3):end-1+dim(3));
    u1 = u(1:end-2,1:end-2,1:end-2);
    u2 = u(1+dim(1):end-2+dim(1)...
          ,1+dim(2):end-2+dim(2)...
          ,1+dim(3):end-2+dim(3));
    u3 = u(1+2*dim(1):end-2+2*dim(1)...
          ,1+2*dim(2):end-2+2*dim(2)...
          ,1+2*dim(3):end-2+2*dim(3));
    v = ((u2-u1).*h2.^2+(u3-u2).*h1.^2)./(h1.*h2.*(h1+h2));
end

