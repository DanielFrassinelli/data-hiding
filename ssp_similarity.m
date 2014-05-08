function [similarity , b] = ssp_similarity(wimg, img, wat, alpha)
%  SSP_SIMILARITY
%    Watermark extraction and similarity measure for multiplicative spread
%    spectrum.
timg = reshape(dct2(img), 1, []);
twimg = reshape(dct2(wimg), 1, []);
wat = reshape(double(wat), 1, []);

timg = sort(timg, 'descend');
twimg = sort(twimg, 'descend');

v = timg(2:1025);
w = twimg(2:1025);

extwat = (w - v) / alpha ./ v;

% ensure that similarity with the very same watermarked image is infinity
similarity = inf;
similarity(any(extwat)) = extwat * wat.' / norm(extwat);

b = wpsnr(wimg , img);

end
