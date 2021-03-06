function threshold = computeThreshold(vmark)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mark_size = 1024;

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

end

function mark = transformMark(mark) % test
    for j = 1:size(mark,2) 
        if (mark(j) == 0)
            mark(j) = -1;
        end
    end
end

