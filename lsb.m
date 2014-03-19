function [watermarked] = lsb(img, watermark, bits)

wx = size(img, 1);
wy = size(img, 2);

watermarked = img;
for i = 1:wx
    for j = 1:wy
        w = bitshift(watermark(i, j), bits - 8, 'uint8');
        p = img(i, j);
        watermarked(i, j) = p - rem(p, 2^bits) + w;
    end
end
