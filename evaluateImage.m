function evaluateImage(original,marked)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

imshow(original);
figure(); 
imshow(marked);

disp('---------- WPSNR ------------');
fprintf('--> %5.2f DB\n' , wpsnr(original,marked));
end

