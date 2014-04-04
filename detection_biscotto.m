function [a, b] = detection_biscotto(imgfile, watfile, attimg)
%  The detection file:
%  - IS UNIQUE for all images,
%  - it must not call any other external function,
%  - must not open any pop-up windows,
%  - must have the following format:
%
% [a, b] = detection_biscotto(imgfile, watfile, attimg)
%
% where
%  - `a` = 0 if the mark is not present or not recognized
%  - `a` = 1 if the mark is recognized in the image under test
%  - `b` = WPSNR between watermarked image and watermarked and attacked image
%  - `attimg` is the string corresponding to the original imagefile
%  - `watfile` is the string corresponding to the watermarked imagefile
%  - `attimg` is the string corresponding to the attacked (under test) imagefile
%
%  The images must be read inside the detection function.
%
% The detection threshold is computed as follows: test your watermarked image
% by presenting 1000 random marks (including the given one) to the detection
% system. Consider the second highest correlation peak value p2. The
% threshold is obtained by: p2+10%(p2). It is unique for all your images and
% it must be computed in advance. It can not be computed during the competion
% in the detection phase.
img = imread(imgfile);
watermark = imread(watfile) / 255.;
attacked = imread(attimg);

watermarked = embedding_biscotto(img, watermark);

a = ssp_similarity(attacked, img, watermark, 0.3) > 6;
b = wpsnr(watermarked / 255., attacked / 255.);

end