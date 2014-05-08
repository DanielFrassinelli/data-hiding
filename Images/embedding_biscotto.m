function [watermarked] = embedding_biscotto(image)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% image : 512x512 unit8
% mark  : 32x32   0,1 double
% alpha : positive double

%% --------- Constants -------- %

alpha = 0.15;    % alpha value (more or less random)

our_seed = 828;

block_size  = 16; % 16 

imsize      = 512; % 512

mark_size   = 1024; %1024 = 32*32

similar_size = 8; % number of similar block we want to use

dct_list_size = imsize*imsize/(block_size * block_size); % number of blocks

%% --------- Preprocessing -------- %

rng(our_seed);   % init prng
mark = round(rand(32,32));
mark = reshape(mark, 1, mark_size); % mark
mark = transformMark(mark);

im_dct = blkproc(double(image),[block_size block_size], 'dct2');

dct_list = zeros(dct_list_size ,block_size, block_size);

k = 1; % from block DCT to list DCT
for i = 1:block_size:imsize
    for j = 1:block_size:imsize
        dct_list(k,:,:) = im_dct(i:(i+block_size-1),j:(j+block_size-1));
        k = k + 1;
    end
end

mydata.param     = 0; % embedding param
mydata.row       = 0; % row of the used freq
mydata.col       = 0; % column
mydata.mean      = 0; % mean of the block
mydata.var       = 0; % variance of the block
mydata.sInd      = zeros(1 , similar_size); 

block = repmat(struct(mydata), 1, dct_list_size);

for i = 1:dct_list_size
    block(i).mean = mean(mean(dct_list(i,:,:)));
    block(i).var  = var(var(dct_list(i,:,:)));   
end

for i = 1:dct_list_size   % calcolo i blocchi piu simili per ogni blocco 
    mi = block(i).mean;
    vi = block(i).var;
    a = [block.var];
    term1 = log(1/4 .* ((a ./ vi) + (vi ./ a + 2)));
    term2 = 1/4 * (((mi - [block.mean]).^2)./(vi - a));
    temp  = (1/4 .* term1) + (1/4 .* term2); 
    temp(i) = inf;
    [~ , ind]  = sort(temp , 'ascend');
    block(i).sInd(:) = ind(1:similar_size);
end

temp = zeros(1 , similar_size);
z    = zeros(1 , 3*3);
for i = 1:dct_list_size
    p = 1;
    for r = 1:3
        for c = 1:3
            for k = 1:similar_size
                temp(k) = dct_list(block(i).sInd(k) , r , c);
            end            
            z(p) = var(temp);
            if (c == 1 && r == 1) 
                z(p) = inf;
            end
            p = p+1;
        end
    end
    [~ , ind] = min(z);
    [a,b] = ind2sub(3,ind);
    block(i).row = a;
    block(i).col = b;
    block(i).param = computeBlockParameter(dct_list(i,:,:), block(i).row, block(i).col); 
end

%% ----------- Embedd the watermark ----------------------------- %

for i = 1:mark_size 
    r = block(i).row;
    c = block(i).col;
    dct_list(i,r,c) = (dct_list(i,r,c)) + (block(i).param * alpha * mark(i));
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

fprintf('wpsnr : %f \n' , WPSNR(image , watermarked));
fprintf('salvare l^immagine come BMP\n');
end

%% ----------------------------- Utils ------------------------------------------ %%

function param = computeBlockParameter(block, x , y)
    temp_v = reshape(block , size(block,2) , size(block,3));    
    param = sqrt(norm(temp_v)) + 5 * x * y;
end

function mark = transformMark(mark) %still testing
    for j = 1:size(mark,2) 
        if (mark(j) == 0)
            mark(j) = -1;
        end
    end
end