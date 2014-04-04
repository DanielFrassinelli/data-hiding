function [similarity] = ssp_similarity(wimg, img, wat, alpha)
%  SSP_SIMILARITY
%    Watermark extraction and similarity measure for multiplicative spread
%    spectrum.
timg = reshape(dct2(img), 1, []);
twimg = reshape(dct2(wimg), 1, []);

timg = sort(timg, 'descend');
twimg = sort(twimg, 'descend');

v = timg(2:1001);
w = twimg(2:1001);

extwat = (w - v) / alpha ./ v;
similarity = extwat * wat.' / norm(extwat);

end
