function compressed = applyJPEG(image , rate)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

imwrite(image , 'temp.jpg', 'jpg' , 'Quality' , rate);
compressed = imread('temp.jpg');

end

