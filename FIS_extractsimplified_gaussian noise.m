clc;
clear all;
close all;
% save start time
start_time=cputime;

% read in the cover object

%file_name='Lena256.bmp';
file_name1='Cameraman256.bmp';
I=imread(file_name1);
M=0;V=0.05;
file_name = imnoise(I,'gaussian',M,V);

cover_object1=file_name;
cover_object=im2double(cover_object1);
% determine size of cover image

Mc=size(cover_object,1);	        %Height
Nc=size(cover_object,2);	        %Width

% blocksize=8
blocksize=8; 

% determine maximum message size based on cover object, and blocksize
max_message=Mc*Nc/(blocksize^2);

% process the image in blocks
x=1;
y=1;
sum_dc=0;
  for kk=1:max_message
       
        % transform block 
        image_block=dct2(cover_object(y:y+blocksize-1,x:x+blocksize-1));
        sum_dc=sum_dc + image_block(1,1);
    if (x+blocksize) >= Nc
        x=1;
        y=y+blocksize;
    else
        x=x+blocksize;
    end
        
        end
%claculate mean of all dc cofficients
mean_dc=sum_dc/max_message;

        % process the image in blocks
x=1;
y=1;
sum_dc=0;
  for kk=1:max_message
       
        % transform block 
        image_block=dct2(cover_object(y:y+blocksize-1,x:x+blocksize-1));
        % calculate luminance senstivity
        luminance(1,kk)=(image_block(1,1)/mean_dc)/2;     %------- range(0-2)
        % calculate average contrast
        [te,varn]=statxture(image_block);
        %average_contrast(1,kk)=te(2);
        % calculate threshold
        threshold(1,kk)=graythresh(image_block);        %-------range(0-1)
        % calculate normalized variance value
        variance(1,kk)=varn;                            %------- range(0-1)
        % move on to next block. At and of row move to next row
        if (x+blocksize) >= Nc
            x=1;
            y=y+blocksize;
        else
            x=x+blocksize;
        end
        
  end
%FIS system for generating watermark invoked & 
%inputs given to it for each block 
dct_fuz=readfis('watermark_2');  
for kk=1:max_message
    w(1,kk)=evalfis([ luminance(1,kk) threshold(1,kk) variance(1,kk)],dct_fuz);
end

%Whole 256x256 image's dct taken
dct_cover=dct2(cover_object);
val=dct_cover(1,1);
dct_cover(1,1)=min(min(dct_cover)); %for 2d matrix 1st column wise then row 
watermark=randn(1,max_message);

%sort the DCT of cover image so as to obtain the low frequency components
[svals,idx] = sort(dct_cover(:),'descend'); % sort to vector
lvals=svals;
k=0.05;
%idx returns the corresponding indices of all unsorted elements
%useful in preserving the orignal location of all elements

%Take top kk DCT cofficients and spread the watermark (noise)
for i=1:max_message
    svals(i)=svals(i) + (k * watermark(1,i) * w(1,i));
end



% store position in matrix of top 500000 DCT cofficients
for i=1:max_message
[II,JJ] = ind2sub([Mc,Nc],idx(i)); % position in the matrix
row(i)=II;
col(i)=JJ;
end

% transform the sorted vector again into matrix form
for i=1:(Mc*Nc)
        [II,JJ] = ind2sub([Mc,Nc],idx(i));
        dct_watermark(II,JJ)=svals(i);
end


dct_watermark(1,1)=val;
%Take the inverse DCT
watermarked_image=idct2(dct_watermark);

% convert to uint8 and write the watermarked image out to a file
watermarked_image_int=im2uint8(watermarked_image);
imwrite(watermarked_image_int,'dct_fuzzy.bmp','bmp');

% display processing time
elapsed_time=cputime-start_time,

imshow(watermarked_image,[])

i=imread('Cameraman256.bmp');
j=imread('dct_fuzzy.bmp');
psnr1(i,j);
dlmwrite('dct_fuzzywatermark.txt',watermark);
dlmwrite('dct_fuzzyrow.txt',row);
dlmwrite('dct_fuzzycol.txt',col);


%DECODING THE IMAGE AND COMPARING BOTH PIXEL VALUES by cox method

for i=1:max_message
diff(i)=svals(i)-lvals(i);
end

%Extracted Watermark = (DCT_Low_Frequency_cofficients_Of_Signed_image -
%DCT_Low_Frequency_cofficients_Of_Original image) / (k * weighting factor)

for i=1:max_message
    watermark_n(1,i)=diff(i)/(k*w(1,i));
end
watermark_n(1,1)=0; %dc component is unchanged

SIM=sum(watermark.*watermark_n)/sqrt(sum(watermark.*watermark_n)); 

SIM

