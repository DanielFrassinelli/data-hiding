function [watermarked] = spread_spectrum(img, watermark, alpha)
% SPREAD_SPECTRUM Multiplicative Spread Spectrum watermarking.
%
%  watermarked = SPREAD_SPECTRUM(img, watermark, alpha)
%      Applies a 32x32 matrix of pseudorandom bits `watermark` to the
%      image `img`, of 512x512 pixels, with intensity `alpha`.
watermark = reshape(watermark, 1, []);
timg = reshape(dct2(img), 1, []);
stimg = sign(timg);
atimg = abs(timg);
[sorted, indexes] = sort(atimg, 'descend');

awatermarked = atimg;
for j = 1:1024
    i = indexes(j+1);
    awatermarked(i) = atimg(i) * (1 + alpha * watermark(j));
end

watermarked = awatermarked .* stimg;
watermarked = reshape(watermarked, size(img, 1), size(img, 2));
watermarked = uint8(idct2(watermarked));

end
