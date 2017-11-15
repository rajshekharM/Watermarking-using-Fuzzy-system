function PSNR=psnr1(i,j)
%j=imread('dct_fuzzy.bmp');
%i=imread('Lena256.bmp');

%I=im2double(i);
%J=im2double(j);
error =imsubtract(i,j);
MSE= mean(error(:).^2)
PSNR = 10 * log10(255^2 / MSE)
%PSNR = round(PSNR*100)/100

 
end