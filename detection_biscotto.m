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

alpha = 0.15; % test value

our_seed = 828; % seed

block_size  = 16; % 16 - 8

imsize      = 512; % 512

mark_size   = 1024; %1024 = 32*32

similar_size = 8; % number of similar block we want to use

dct_list_size = imsize*imsize/(block_size * block_size); % number of blocks

threshold = 5.0; % now is dynamic --> we will fix it during the competition

%% --------- Preprocessing-------- %

original = imread(imgfile);         % original image
marked =   imread(attimg);            % filtered image
marked_original  = imread(watfile); % watermarked image

%original = imgfile;         % original image
%marked =   attimg;            % filtered image
%marked_original  = watfile; % watermarked image

rng(our_seed);   % init prng
mark = round(rand(32,32));
mark = reshape(mark, 1, mark_size); % mark
mark = transformMark(mark);

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

computed_param = zeros(dct_list_size);                   % embedding parameter
computed_index = zeros(dct_list_size , 2);               % index of the used frequency
%similarBlockIndex = zeros(dct_list_size , similar_size); % list of x similar block ordered by similarity

for i = 1:dct_list_size   
    temp = reshape(original_dct_list(i,:,:), 1 , block_size*block_size);
    [~ , ind] = sort(abs(temp) , 'descend');
    r = mod(ind(2) - 1,16) + 1;
    c = ceil(ind(2)/16);
    
    computed_param(i) = computeBlockParameter(original_dct_list(i,:,:), r, c);
    computed_index(i,:) = [r,c];
end
%{
values = zeros(dct_list_size , 2);
    
for i = 1:dct_list_size
    values(i, :) = [mean(mean(marked_original_dct_list(i,:,:))) , var(var(marked_original_dct_list(i,:,:)))];
end    

for i = 1:dct_list_size
    temp = zeros(1 ,dct_list_size);
    m = values(i , 1);
    
    for j = 1:dct_list_size
        temp(j) = abs(values(j , 1) - m) * abs(marked_original_dct_list(i,1,2) - marked_original_dct_list(j,1,2));
    end
    
    [~ , ind] = sort(temp , 'ascend');
    similarBlockIndex(i,:) = ind(2:1+similar_size);

    if (i <= 5)
        fprintf('Block : %d ---------------------------------------------------- \n' , i);
        
        base_freq = abs(marked_original_dct_list(i,1,2));
        filt_freq = abs(marked_dct_list(i,1,2));

        rapp_freq = (filt_freq - base_freq)/ base_freq; % variazione frequenza
        
        fprintf('[%4d\t %9.4f\t %9.4f\t %f]\n' , i , marked_original_dct_list(i,1,2) , marked_dct_list(i,1,2) , rapp_freq);
        
        result = 0;
        %t1  = (-((2.5*base_freq)^2))/log(0.1);
        t2  = (-(20)^2)/log(0.5);
        
        for j = 1:similar_size       
            index = similarBlockIndex(i,j);
                        
            sbase_freq = abs(marked_original_dct_list(index,1,2));
            sfilt_freq = abs(marked_dct_list(index,1,2));
        
            s1 = abs(sbase_freq / base_freq);
            s2 = 1 - exp((-(sbase_freq)^2)/t2);
            
            if (s1 > 1)
                s1 = 1/s1;
            end
            
            f = (sfilt_freq - sbase_freq) / sbase_freq; 
            f = f * s1 * s2;
            result = result + f;
        
            fprintf('[%4d\t %9.4f\t %9.4f\t %f] {sim1:%f} {sim2:%f}\n', index , marked_original_dct_list(index,1,2) , marked_dct_list(index,1,2) , f, s1 , s2);
        
        end        
        fprintf('average --> %f\n' , result / similar_size);
    end
end
%}

%% ----------- Exctract the watermark ----------------------------- %

k = 1;
vmark = mark;
for j = 1:mark_size    
    r = computed_index(k,1);
    c = computed_index(k,2);   
    freq = original_dct_list(k,r,c); 
    param = computed_param(k);
    vmark(j)   = (marked_dct_list(k,r,c) - (freq))/(alpha * param); %remove watermark
    k = k + 1;
end

%% --------- Postprocessing -------- %

sim = (vmark(1,:) * mark(1,:)')/norm(vmark); % compute similarity
newsim = zeros(1,1000);
rng('shuffle');

for i = 1:1000
    newmark = round(rand(32,32)); %similary test
    newmark = reshape(newmark, 1, mark_size); % mark
    newmark = transformMark(newmark);
    newsim(i) = (vmark(1,:) * newmark(1,:)')/norm(vmark);
end

ss = sort(newsim , 'descend');
threshold = ss(1) + ss(1)/10.0;
fprintf('\nSimilarity   : %f    Threshold   : %f  \n' , sim , threshold);

if (sim >= threshold)
    a = 1;
else
    a = 0;
end

b = wpsnr(marked_original , marked);

end

%% ----------------------------- Utils ------------------------------------------ %%

function param = computeBlockParameter(block, x , y)
temp_v = reshape(block , size(block,2) , size(block,3));
param = sqrt(norm(temp_v));
end

function dist = computeDistBetweenBlock(a , b , r , c)
a = reshape(a ,16,16);
b = reshape(b ,16,16);
d = a - b;
d(r,c) = d(r,c)*d(r,c);
dist = norm(d);
end

function mark = transformMark(mark) % test
    p = 1.0;
    for j = 1:size(mark,2) 
        if (mark(j) == 0)
            mark(j) = -1;
        end
        if (mod(j,2) == 0)
            mark(j) = mark(j) * p;
            if (p >= 2.0)
                p = 1.0;
            end
            p = p + 0.2;
        end
    end
end

