direc = fullfile('data', 'superdarn');
pngs = {dir(fullfile(direc, '*.png')).name};
wdth = 400;
hght = 300;
padw = 16;
padh = 32;
im_full = uint8(255 * ones((hght + padh) * 3, wdth * 2 + padw, 3));
crop_h = [370, 380, 300, 375, 270, 340];
crop_w = [175, 195, 85, 190, 70, 110];
mlen = 16;
mwdt = 3;
mclr = [100, 0, 255];
moff = 60;
bh = 18;
bw = 18;
bh0s = [270, 270, 95, 270, 125, 50];
bw0s = [0, 0, 15, 0, 30, 0];
toff = -1;
ftsz = 12;
ftsz = round(ftsz * (size(im_full, 2) / 96) / 6.5);
lh0 = 5;
lw0 = 5;

for i = 1:6
    png = pngs{i};
    time = datetime(png(1:end-4), 'InputFormat', 'uuuuMMdd_HHmmss');
    time.Format = 'MMM dd, HH:mm';
    im = imread(fullfile(direc, png));
    
    % legend = im(647:682, 675:721, :);
    legend1 = im(647:682-29, 675:721, :);
    legend2 = im(647+20:682, 675:721, :);
    legend = [legend1; legend2];

    ch0 = crop_h(i) + 1;
    cw0 = crop_w(i) + 1;
    im_crop = im(ch0:ch0+hght-1, cw0:cw0+wdth-1, :);

    mh0 = round(hght/2 - mlen/2) + 1 + moff;
    mw0 = round(wdth/2 - mlen/2) + 1;
    mh1 = round(hght/2 - mwdt/2) + 1 + moff;
    mw1 = round(wdth/2 - mwdt/2) + 1;
    mline = permute(repmat(mclr, [mwdt, 1, mlen]), [1, 3, 2]);
    im_crop(mh1-mlen/2:mh1-mlen/2+mwdt-1, mw0:mw0+mlen-1, :) = mline;
    im_crop(mh1+mlen/2:mh1+mlen/2+mwdt-1, mw0:mw0+mlen-1, :) = mline;
    im_crop(mh0:mh0+mlen-1, mw1-mlen/2:mw1-mlen/2+mwdt-1, :) = permute(mline, [2, 1, 3]);
    im_crop(mh0:mh0+mlen-1, mw1+mlen/2:mw1+mlen/2+mwdt-1, :) = permute(mline, [2, 1, 3]);
    
    bh0 = bh0s(i) + 1;
    bw0 = bw0s(i) + 1;
    im_crop(bh0:bh0+bh-1, bw0:bw0+bw-1, :) = ones(bh, bw, 3) * 255;
    
    im_crop(hght-size(legend, 1)+1-lh0:hght-lh0, 1+lw0:size(legend, 2)+lw0, :) = legend;

    im_crop(1, :, :) = ones(1, wdth, 3) * 128;
    im_crop(end, :, :) = ones(1, wdth, 3) * 128;
    im_crop(:, 1, :) = ones(hght, 1, 3) * 128;
    im_crop(:, end, :) = ones(hght, 1, 3) * 128;

    im_crop = insertText(im_crop, [0, 0], char(64+i), ...
        'AnchorPoint', 'LeftTop', 'BoxOpacity', 0, 'FontSize', ftsz);

    h0 = ceil(i/2 - 1) * (hght + padh) + 1 + padh;
    w0 = mod(i-1, 2) * (wdth + padw) + 1;
    im_full(h0:h0+hght-1, w0:w0+wdth-1, :) = im_crop;
    
    th = round(padh / 2) + ceil(i/2 - 1) * (hght + padh) + 1 - toff;
    tw = round(wdth / 2) + mod(i-1, 2) * (wdth + padw) + 1;
    labl = sprintf('%s UT', time);
    im_full = insertText(im_full, [tw, th], labl, ...
        'AnchorPoint', 'Center', 'TextBoxColor', 'white', 'FontSize', ftsz);
end

% imshow(im_full)
imwrite(im_full, 'plots/00_superdarn_summary.png', 'png')
% close all