clc
clear
close all

dirinput = 'testimages/blastocysts';
diroutput = 'testoutputs/blastocysts';

if ~exist(diroutput,'dir')
  mkdir(diroutput);
end

filesinput = dir(dirinput);
fileindinput = find(~[filesinput.isdir]);

for i = 1:numel(fileindinput)
  filenameinput = fullfile(dirinput,filesinput(fileindinput(i)).name);
  fprintf('\nSegmenting %s: Image %d of %d',...
    filesinput(fileindinput(i)).name,i,numel(fileindinput));
  inimage = imread(filenameinput);
  imshow(inimage);title('Original Image');pause(1)

  fileprefix = split(filesinput(fileindinput(i)).name,'.');

  inimage = remove_borders(inimage);
  inimage = remove_scale(inimage);
%   imshow(inimage);pause

  if size(inimage,3) == 3
    inimage = rgb2gray(inimage);
  end

  threshimage = thresh_sobel(inimage);
  squareimage = square_bound(inimage, threshimage);
  imshow(squareimage);title('Cropped Image'); pause(2)

  outname = fullfile(diroutput,filesinput(fileindinput(i)).name);
  imwrite(squareimage,outname);
end

fprintf('\nEnd of Program');

%%                              Removal Functions
% ======================================================================= %
                 %%%  Remove black borders in images %%%
function out_img = remove_borders(in_img)
  out_img = in_img(4:end-2,4:end-7,:);
end

                      %%% Remove green scale bar %%%
function out_img = remove_scale(in_img)
  % Scan through corner points
  for i = round(size(in_img,1)*7/8):size(in_img,1)-5
    for j = 400:size(in_img,2)-2
      % Perform analysis on green pixels
      if in_img(i,j,2) == 255 && in_img(i,j,1) == 0 &&...
          in_img(i,j,3) == 0
        for k = 1:3
          in_img(i,j,k) = mean_window(in_img,i,j,k);
        end
      end
    end
  end
  out_img = in_img;
end

function out_pixel = mean_window(in_img,r,c,channel)
  window = in_img(r-2:r+2,c-2:c+2,channel);
  window_filt = window(window~=0);
  out_pixel = mean(window_filt);
end

%%          Thresholding using Sobel and Morphological Operations
% ======================================================================= %
function out_img = thresh_sobel(in_img)
  [~,threshold] = edge(in_img,'sobel');
  fudgeFactor = 0.5;

  BWs = edge(in_img,'sobel',threshold * fudgeFactor);
  se90 = strel('line',3,90);
  se0 = strel('line',3,0);
  BWsdil = imdilate(BWs,[se90 se0]);

  BWdfill = imfill(BWsdil,'holes');
  BWCC = bwareafilt(BWdfill,1);
  BWCONV = bwconvhull(BWCC);
  
  DIFcont = logical(BWCONV - BWCC);
  sesq = strel('disk',20);
  DIFcont = imclose(DIFcont,sesq);

  BWdif = BWsdil & imcomplement(DIFcont);
  BWCC = bwareafilt(BWdif,1);
  BWCC = imfill(BWCC,'holes');

  sesq = strel('diamond',5);
  BWCC = imclose(BWCC,sesq);
  out_img = bwconvhull(BWCC);
end 

%%           Crop Image to Smallest Bounding Square using Mask
% ======================================================================= %
function out_img = square_bound(in_img,thresh_img)
  [m,n] = size(in_img);
  super_thresh = padarray(thresh_img,[m n]);
  super_img = padarray(in_img,[m n],'symmetric');

  serad = round(0.20*sqrt(bwarea(thresh_img)/pi));
  sedil = strel('disk',serad);
  super_thresh = imdilate(super_thresh,sedil);

  stats = regionprops(super_thresh,'Centroid','Orientation',...
    'MajorAxisLength','MinorAxisLength');
  Cx = stats.Centroid(1);
  Cy = stats.Centroid(2);
  Rx = stats.MajorAxisLength/2;
  Ry = stats.MinorAxisLength/2;
  if stats.Orientation < 0
      theta = pi - (stats.Orientation + 180)*pi/180;
  else
      theta = pi - stats.Orientation*pi/180;
  end

  ellip = zeros(size(super_thresh));
  ellip = imbinarize(ellip);
  t = linspace(0,2*pi,500);
  x = round(Rx*cos(t)*cos(theta) - Ry*sin(t)*sin(theta) + Cx);
  y = round(Rx*cos(t)*sin(theta) + Ry*sin(t)*cos(theta) + Cy);
  ind = sub2ind(size(ellip),y,x);
  ellip(ind) = 1;
  rad = 5;
  while 1
    SE = strel('disk',rad);
    ellip = imdilate(ellip,SE);
    ellip = imfill(ellip,'holes');
    stats = regionprops(ellip,'Area');
    if length(stats) == 1 
      break
    end
    rad = rad + 5;
  end

  stats = regionprops(ellip,'BoundingBox');
  left = ceil(stats.BoundingBox(1));
  top = ceil(stats.BoundingBox(2));
  width = stats.BoundingBox(3);
  height = stats.BoundingBox(4);

  if mod((top+height),2) == 1
    height = height-1;
  end
  if mod((left+width),2) == 1
    width = width-1;
  end
  if height ~= width
    if height == max([height,n])
      dif = height - width;
      left = left-round(dif/2);
      width = width+dif;
    else
      dif = width - height;
      top = top-round(dif/2);
      height = height+dif;
    end
  end

  out_img = super_img(top:top+height,left:left+width);
end

