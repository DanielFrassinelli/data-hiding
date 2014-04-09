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

% --> al momento imgfile è la matrice dell'immagine originale e watfile è la matrice
% del watermark ... cambiero più avanti, attacked è la matrice
% dell'immagine col watermark  filtri

%img = imread(imgfile);
%watermark = imread(watfile) / 255.;
%attacked = imread(attimg);

%% --------- Constants -------- %

alpha = 0.5; % alpha

block_size  = 16; % 16 - 8

imsize      = 512; % 512

mark_size   = 1024; %1024 = 32*32

dct_list_size = imsize*imsize/(block_size * block_size); % number of blocks

threshold = 5.0;

%% --------- Preprocessing -------- %

mark = watfile;
original = imgfile;
marked = attimg;

mark = reshape(mark , 1 ,mark_size);

for j = 1:mark_size
    if (mark(j) == 0)
        mark(j) = -1;
    end
end

marked_original = embedding_biscotto(original , mark); % here we just need to load the image 

original_dct                = blkproc(double(original),         [block_size block_size], 'dct2');
marked_dct                  = blkproc(double(marked),           [block_size block_size], 'dct2');
marked_original_dct         = blkproc(double(marked_original),  [block_size block_size], 'dct2');

original_dct_list        =  zeros(dct_list_size ,block_size,block_size);
marked_dct_list          =  zeros(dct_list_size ,block_size,block_size);
marked_original_dct_list =  zeros(dct_list_size ,block_size,block_size);

k = 1; % from block DCT to list DCT
for i = 1:block_size:imsize
    for j = 1:block_size:imsize
        original_dct_list(k,:,:) = original_dct(i:(i+block_size-1),j:(j+block_size-1));
        marked_dct_list(k,:,:)  = marked_dct(i:(i+block_size-1),j:(j+block_size-1));
        marked_original_dct_list(k,:,:)  = marked_original_dct(i:(i+block_size-1),j:(j+block_size-1));
        k = k + 1;
    end
end

%% ----------- Decide which blocks we want to use ------------------------- %

computed_param = zeros(dct_list_size);
computed_index = zeros(dct_list_size , 2);

for i = 1:dct_list_size   
    temp = reshape(original_dct_list(i,:,:), 1 , block_size*block_size);
    [~ , ind] = sort(abs(temp) , 'descend');
    r = mod(ind(2) - 1,16) + 1;
    c = ceil(ind(2)/16);
    
    computed_param(i) = computeBlockParameter(original_dct_list(i,:,:), r, c);
    computed_index(i,:) = [r,c];

end

%% ----------- Exctract the watermark ----------------------------- %

k = 1;
vmark = mark;
for j = 1:mark_size    
    r = computed_index(k,1);
    c = computed_index(k,2);
    temp = marked_dct_list(k,r,c) - original_dct_list(k,r,c);
    vmark(j)   = temp/(alpha * computed_param(k)); %remove watermark
    k = k + 1;
end

%% --------- Postprocessing -------- %

sim = (vmark(1,:) * mark(1,:)')/sqrt(vmark(1,:) * vmark(1,:)'); % compute similarity

fprintf('\nSimilarity   : %f\n' , sim);
if (sim >= threshold)
    a = 1;
else
    a = 0;
end

b = wpsnr(marked_original , marked);

end

%% ----------------------------- Utils ------------------------------------------ %%

function param = computeBlockParameter(block, x , y)

temp_v = reshape(block , 1 , size(block , 2) * size(block , 3));
param = sqrt(norm(temp_v));

end
