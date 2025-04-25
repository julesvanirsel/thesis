comparisons = [...
    ["swop_20230210_35487_09_SD_AM_AC"; "swop_20230210_35487_09_PF_AM_AC"], ...
    ["swop_20230210_35487_09_SD_AM_AC"; "swop_20230210_35487_09_SD_UM_AC"], ...
    ["swop_20230304_27012_09_SD_AM_xC"; "swop_20230304_27012_09_NB_AM_xC"], ...
    ["swop_20230314_24547_09_SD_AM_AC"; "swop_20230314_24547_09_PF_AM_AC"], ...
    ];

for comp = comparisons
    clear('sim')
    sim.A = char(comp(1));
    sim.B = char(comp(2));
    fprintf('%s v %s\n', sim.A, sim.B)

    direc_root = getenv('GEMINI_SIM_ROOT');
    resf = 4;
    ftsz = 14 * resf;
    letter_pos = [[0, 0]; [250, 0]; [0, 220]; [250, 220]] * resf;

    first = true;
    for s = 'AB'
        direc = fullfile(direc_root, sim.(s), 'plots3d');

        filename_iso = dir(fullfile(direc, '*ISO_P_000.png')).name;
        filename_side = dir(fullfile(direc, '*SID_P_000.png')).name;
        filename_top = dir(fullfile(direc, '*TOP_P_000.png')).name;

        im_iso = imread(fullfile(direc, filename_iso));
        im_sid = imread(fullfile(direc, filename_side));
        im_top = imread(fullfile(direc, filename_top));

        im_top = imshift(im_top, 0, round(4.5*resf));

        im_sid = im_sid(1:end-34*resf, 1+10*resf:end-0*resf, :);
        im_top = im_top(1:end-34*resf, 1+20*resf:end-0*resf, :);

        im_clb_j = im_iso(1:250*resf, end-75*resf:end, :);
        im_clb_n = im_iso(end-250*resf+1:end, end-75*resf:end, :);
        im_clb_j = imshift(im_clb_j, 0, -5*resf);
        im_clb_n = imshift(im_clb_n, 0, -9*resf);

        im_lbl_in = im_iso(19*resf:32*resf, 22*resf:176*resf, :);
        im_lbl_out = im_iso(round((19+15.5)*resf):round((32+15.5)*resf), ...
            22*resf:176*resf, :);

        im_top(1+8*resf:size(im_lbl_in, 1)+8*resf, ...
            1+44*resf:size(im_lbl_in, 2)+44*resf, :) = im_lbl_in;
        im_top(1+8*resf:size(im_lbl_out, 1)+8*resf, ...
            1+198*resf:size(im_lbl_out, 2)+198*resf, :) = im_lbl_out;

        cntr_crop = 34 * resf;
        if first
            im = uint8(ones(size(im_sid, 1) * 2, ...
                size(im_sid, 2) + size(im_top, 2) + size(im_clb_j, 2), 3)) * 255;
            im_sid = im_sid(1:end-cntr_crop, :, :);
            im_top = im_top(1:end-cntr_crop, :, :);
            im(1:size(im_sid, 1), 1:size(im_sid, 2), :) = im_sid;
            im(1:size(im_sid, 1), ...
                size(im_sid, 2)+1:size(im_sid, 2)+size(im_top, 2), :) = im_top;
            first = false;
        else
            im(size(im_sid, 1)+1-cntr_crop:end-cntr_crop, 1:size(im_sid, 2), :) = im_sid;
            im(size(im_sid, 1)+1-cntr_crop:end-cntr_crop, ...
                size(im_sid, 2)+1:size(im_sid, 2)+size(im_top, 2), :) = im_top;
        end
    end

    im = im(1:end-cntr_crop, :, :);
    im(1:size(im_clb_j, 1), end-size(im_clb_j, 2)+1:end, :) = im_clb_j;
    im(end-size(im_clb_n, 1)+1:end, end-size(im_clb_n, 2)+1:end, :) = im_clb_n;

    for i = 1:length(letter_pos)
        im = insertText(im, letter_pos(i, :), char(64 + i), ...
            'AnchorPoint', 'LeftTop', 'BoxOpacity', 0, 'FontSize', ftsz);
    end

    hght = ceil(size(im, 1) * 2 / 96 / resf) *96 * resf / 2;
    wdth = ceil(size(im, 2) * 2 / 96 / resf) *96 * resf / 2;
    im_new = uint8(ones(hght, wdth, 3)) * 255;
    im_new(1:size(im, 1), 1:size(im, 2), :) = im;

    filename = sprintf('plots/00_gemini_compare_%s_v_%s.png', sim.A(6:end), sim.B(24:end));
    imwrite(im_new, filename, 'png')

end

function im_shift = imshift(im, right, down)
im_shift = circshift(circshift(im, right, 2), down, 1);
end
