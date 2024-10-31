x=linspace(-500, 500, 400); % km
y=linspace(-200, 200, 600); % km
[X, Y] = ndgrid(x, y);

dx = median(diff(x));
dy = median(diff(y));

s1 = 60; % km
s2 = 15; % km
smin = 15; % km

f1 = 2 * pi / s1;
f2 = 2 * pi / s2;

wx = 2 * smin / dx;
wy = 2 * smin / dy;

close all
data = sin(f1*X).*sin(f1*Y);
data = data + 0*sin(f2*X).*sin(f2*Y) / 10;
data = data + 0*rand(size(X)) / 20;
data_smooth = smoothdata2(data, 'gaussian', {wx, wy});

ar = [range(x), range(y), 300];

figure('Position', [-1800, 50, 800, 800])
pcolor(X, Y, data_smooth ./ data)
pbaspect(ar)
colorbar
zlim([-1,1])
clim([-1,1])
shading flat
