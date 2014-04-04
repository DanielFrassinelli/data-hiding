function [watermarked] = spread_spectrum(img, watermark, alpha)
% SPREAD_SPECTRUM Multiplicative Spread Spectrum watermarking.
%
%  watermarked = SPREAD_SPECTRUM(img, watermark, alpha)
%      Applies a 1000-len vector of pseudorandom bits `watermark` to the
%      image `img`, of at least 1000 pixels, with intensity `alpha`.

timg = reshape(dct2(img), 1, []);
stimg = sign(timg);
atimg = abs(timg);
[sorted, indexes] = sort(atimg, 'descend');

awatermarked = atimg;
for j = 1:1000
    i = indexes(j+1);
    awatermarked(i) = atimg(i) * (1 + alpha * watermark(j));
end

watermarked = awatermarked .* stimg;
watermarked = reshape(watermarked, size(img, 1), size(img, 2));
watermarked = uint8(idct2(watermarked));

end
