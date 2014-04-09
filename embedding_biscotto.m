function [watermarked] = embedding_biscotto(image, mark)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% image : 512x512 unit8
% mark  : 32x32   0,1 double
% alpha : positive double

%% --------- Constants -------- %

alpha = 0.5;    % alpha value, (preso a caso al momento)

block_size  = 16; % 16 

imsize      = 512; % 512

mark_size   = 1024; %1024 = 32*32

dct_list_size = imsize*imsize/(block_size * block_size); % number of blocks

%% --------- Preprocessing -------- %

mark = reshape(mark , 1 ,mark_size);

for j = 1:mark_size
    if (mark(j) == 0)
         mark(j) = -1;
    end
end

im_dct = blkproc(double(image),[block_size block_size], 'dct2');

dct_list = zeros(dct_list_size ,block_size, block_size);

k = 1; % from block DCT to list DCT
for i = 1:block_size:imsize
    for j = 1:block_size:imsize
        dct_list(k,:,:) = im_dct(i:(i+block_size-1),j:(j+block_size-1));
        k = k + 1;
    end
end

%% ----------- Compute parameters we want to use ------------------------- %

computed_param = zeros(dct_list_size);
computed_index = zeros(dct_list_size , 2);

for i = 1:dct_list_size   
    temp = reshape(dct_list(i,:,:), 1 , block_size*block_size);
    [~ , ind] = sort(abs(temp) , 'descend');
    r = mod(ind(2) - 1,16) + 1;
    c = ceil(ind(2)/16);

    computed_param(i) = computeBlockParameter(dct_list(i,:,:), r, c);
    computed_index(i,:) = [r,c];

end

%% ----------- Embedd the watermark ----------------------------- %

k = 1;
for j = 1:mark_size 
    r = computed_index(k,1);
    c = computed_index(k,2);
   
    dct_list(k,r,c) = (dct_list(k,r,c)) + (computed_param(k) * alpha * mark(j));
    k = k + 1;   
end

%% --------- Postprocessing -------- %

watermarked = zeros(imsize,imsize);

k = 1; % from DCT list to DCT block
for i = 1:block_size:imsize
    for j = 1:block_size:imsize
        watermarked(i:(i+block_size-1),j:(j+block_size-1)) = dct_list(k,:,:);
        k = k + 1;
    end
end

watermarked =  uint8(blkproc(watermarked,[block_size block_size], 'idct2'));

fprintf('wpsnr : %f \n' , wpsnr(image , watermarked));

end

%% ----------------------------- Utils ------------------------------------------ %%

function param = computeBlockParameter(block, x , y)

temp_v = reshape(block , 1 , size(block , 2) * size(block , 3));
param = sqrt(norm(temp_v));

end
