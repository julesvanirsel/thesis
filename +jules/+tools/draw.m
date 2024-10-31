function draw(image)
arguments
    image (:, :) {isnumeric}
end

chars = [1, 185, 97, 181, 96, 176, 184, 127, 248, 125, 47, 40, 169, ...
    48, 177, 173, 44, 62, 63, 61, 42, 43, 178, 59, 124, 45, 35, 216, 95, ...
    182, 92, 162, 186, 38, 167, 180, 93, 172, 179, 188, 115, 64, 192, 187, ...
    100, 171, 65, 119, 56, 232, 71, 116, 90, 85, 237, 238, 175, 121, 106, ...
    68, 239, 222, 109, 250, 251, 206, 52, 105, 252, 87, 200, 51, 253, 81, ...
    164, 223, 122, 242, 91, 73, 120, 166, 201, 202, 226, 86, 89, 190, 114, ...
    245, 163, 113, 247, 108, 199, 219, 246, 193, 256, 230, 69, 80, 165, 195, ...
    221, 83, 191, 255, 170, 197, 211, 212, 224, 168, 213, 196, 88, 198, 57, ...
    215, 79, 37, 214, 49, 104, 82, 39, 210, 183] - 1;

chars = chars;

image = image';
image_wdth = size(image, 1);
image_hght = size(image, 2);
wdth = 500;
hght = round(image_hght * wdth / image_wdth);
image_norm = 255 * (image - min(image(:))) / (max(image(:)) - min(image(:)));
fimage = griddedInterpolant(image_norm);

for j = 1:hght
    for i = 1:wdth
        pixel = fimage(round(i*image_wdth/wdth), round(j*image_hght/hght));
        if isnan(pixel)
            id = 1;
        else
            [~, id] = min(abs(pixel - chars));
        end
        % fprintf('%.0f',id)
        fprintf('%s', char(chars(id)))
    end
    fprintf('\n')
end
