sims = [...
    "swop_20230210_35487_09_SD_AM_AC", ...
    "swop_20230210_35487_09_PF_AM_AC", ...
    "swop_20230210_35487_09_SD_UM_AC", ...
    "swop_20230304_27012_09_SD_AM_xC", ...
    "swop_20230304_27012_09_NB_AM_xC", ...
    "swop_20230314_24547_09_SD_AM_AC", ...
    "swop_20230314_24547_09_PF_AM_AC", ...
    ];

for s = sims
    sim = char(s);
    fprintf('%s\n', sim)
    save_iso = true;
    direc_root = getenv('GEMINI_SIM_ROOT');
    direc = fullfile(direc_root, sim, 'plots3d');
    resf = 4;
    ftsz = 14 * resf;
    letter_pos = [[-2, -5]; [-2, 460]; [270, 460]] * resf;

    filename_iso = dir(fullfile(direc, '*ISO_P_000.png')).name;
    filename_side = dir(fullfile(direc, '*SID_p_000.png')).name;
    filename_top = dir(fullfile(direc, '*TOP_p_000.png')).name;

    im_iso = imread(fullfile(direc, filename_iso));
    im_sid = imread(fullfile(direc, filename_side));
    im_top = imread(fullfile(direc, filename_top));

    im_top = imshift(im_top, -10*resf, round(4.5*resf));
    im_iso = imshift(im_iso, -5*resf, 20*resf);

    im_iso = im_iso(25*resf:end, :, :);
    im_sid = im_sid(1:end-40*resf, :, :);
    im_top = im_top(1:end-40*resf, :, :);

    im = uint8(ones(size(im_iso, 1) + size(im_sid, 1), ...
        size(im_iso, 2), 3)) * 255;
    im(1:size(im_iso, 1), :, :) = im_iso;
    im(size(im_iso, 1)+1:end, 1:size(im_sid, 2), :) = im_sid;
    im(size(im_iso, 1)+1:end, size(im_sid, 2)+1:end, :) = im_top;

    hght = ceil(size(im, 1) * 2 / 96 / resf) *96 * resf / 2;
    wdth = ceil(size(im, 2) * 2 / 96 / resf) *96 * resf / 2;

    im_new = uint8(ones(hght, wdth, 3) * 255);
    im_new(1:size(im, 1), 1:size(im, 2), :) = im;

    for i = 1:length(letter_pos)
        im_new = insertText(im_new, letter_pos(i, :), char(64 + i), ...
            'AnchorPoint', 'LeftTop', 'BoxOpacity', 0, 'FontSize', ftsz);
    end

    filename = sprintf('plots/00_gemini_output_%s.png', sim(6:end));
    imwrite(im_new, filename, 'png')
    if save_iso
        im_iso = imshift(im_iso, -5*resf, -8*resf);
        imwrite(im_iso, 'plots/00_fluxtube_example.png', 'png')
    end
end

function im_shift = imshift(im, right, down)
im_shift = circshift(circshift(im, right, 2), down, 1);
end
