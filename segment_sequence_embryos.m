clc
clear
close all

addpath('helpercode')

dirinput = 'testimages/sequence_embryos';
diroutput = 'testoutputs/sequence_embryos';

if ~exist(diroutput, 'dir')
    mkdir(diroutput);
end

filesinput = dir(dirinput);
fileindinput = find(~[filesinput.isdir]);

for i = 1:numel(fileindinput)
  file_name = filesinput(fileindinput(i)).name;
  filename_input = fullfile(dirinput,file_name);
  fprintf('\nSegmenting %s: Image %d of %d',...
    filesinput(fileindinput(i)).name,i,numel(fileindinput));
  img = imread(filename_input);
  imshow(img); title('Original Image'); pause(1)

  % Erase Well Number and TimeStamp
  img(470:end,:) = 0;

  %% Get Micro Well Boundary
  img_pixels = img(:); 
  img_pixels(img_pixels <= 10) = [];
  img_pixels(img_pixels >= 250) = [];

  med_pix = median(img_pixels);
  std_val = round(std(double(img_pixels)));

  mask = img;
  mask(img > (med_pix + std_val)) = 0;
  mask(mask < (med_pix - std_val)) = 0;
  mask(mask > 0) = 1;
  mask = logical(mask);

  mask = imcomplement(mask);
  [M,N] = size(img);
  circ_mask = false(M,N);
  [xx,yy] = meshgrid(1:N,1:M);
  if i ~= 2
      rad = 170;
  else
      rad = 150;
  end
  circ_mask = circ_mask | hypot(xx - N/2, yy - M/2) <= rad;
  mask = mask & ~circ_mask;

  mask_1st = bwareafilt(mask,1);
  se_len = 5; 
  se = strel('square',se_len);
  mask_1st = imerode(mask_1st,se);
  mask_1st = bwareafilt(mask_1st,1,'largest',4);
  mask_1st = imdilate(mask_1st,se);
  mask_1st = mask_1st & bwareafilt(mask,1);


  mask_2nd = bwareafilt(mask,2,'largest',4) & ~bwareafilt(mask,1);
  mask_L2 = mask_2nd;
  mask_L2(:,251:end) = 0;
  mask_L2 = bwareafilt(mask_L2,1);
  mask_R2 = mask_2nd;
  mask_R2(:,1:250) = 0;
  mask_R2 = bwareafilt(mask_R2,1);

  [rowL2,colL2] = ind2sub(size(img),find(mask_L2));
  [min_rowL,~] = max(rowL2); 
  [rowR2,colR2] = ind2sub(size(img),find(mask_R2));
  [min_rowR,~] = max(rowR2);

  if abs(min_rowR - min_rowL) > 100
    disp(file_name);
    disp(abs(min_rowR - min_rowL));
    if min_rowR > min_rowL
      mask_left = mask & ~bwareafilt(mask,2,'largest',4);
      mask_left(:,251:end) = 0;
      mask_left_label = bwlabel(mask_left);
      max_area = 0;
      max_ind = 0;
      for k = 1:max(mask_left_label(:))
        left_comp = mask_left_label == k;
        mask_combined = fliplr(left_comp) & mask_2nd;
        new_area = bwarea(mask_combined);
        if new_area > max_area
          max_area = new_area;
          max_ind = k;
        end
      end
      if max_ind > 0
        left_comp = mask_left_label == max_ind;
        mask_2nd_temp = bwareafilt(mask_2nd,1);
        mask_2nd_temp(251:end,:) = 0;
        [x,y] = minObjectDist(left_comp, mask_2nd_temp);
        mask_2nd = mask_2nd | left_comp;
        mask_2nd(sub2ind(size(img),x,y)) = 1;
      end


    else
      mask_right = mask & ~bwareafilt(mask,2,'largest',4);
      mask_right(:,1:250) = 0;
      mask_right_label = bwlabel(mask_right);
      max_area = 0;
      max_ind = 0;
      for k = 1:max(mask_right_label(:))
        right_comp = mask_right_label == k;
        mask_combined = fliplr(right_comp) & mask_2nd;
        new_area = bwarea(mask_combined);
        if new_area > max_area
            max_area = new_area;
            max_ind = k;
        end
      end
      if max_ind > 0
        right_comp = mask_right_label == max_ind;

        mask_2nd_temp = bwareafilt(mask_2nd,1);
        mask_2nd_temp(251:end,:) = 0;
        [x,y] = minObjectDist(right_comp, mask_2nd_temp);
        mask_2nd = mask_2nd | right_comp;
        mask_2nd(sub2ind(size(img),x,y)) = 1;
      end

    end
  end

  mask = mask_1st | mask_2nd;

  %% Get Inner Area Mask
  mask_perim = bwperim(mask);
  mid_pix = find(mask_perim(:,250));
  mid_pix(mid_pix == 1) = [];
  mid_pix(mid_pix == 500) = [];
  if length(mid_pix) <=2
    mask_end = imclose(mask,strel('square',13));
  else
    mask_L2 = mask_2nd;
    mask_L2(:,251:end) = 0;
    mask_L2 = bwareafilt(mask_L2,1);
    mask_R2 = mask_2nd;
    mask_R2(:,1:250) = 0;
    mask_R2 = bwareafilt(mask_R2,1);

    mask_L1 = mask_1st;
    mask_L1(:,251:end) = 0;
    mask_R1 = mask_1st;
    mask_R1(:,1:250) = 0;

    [rowL2,colL2] = ind2sub(size(img),find(mask_L2));
    [min_rowL,min_indL] = max(rowL2);
    min_colL = colL2(min_indL);
    mask_L1(1:min_rowL,:) = 0;

    [x,y,min_distL] = minPointObjectDist(min_rowL,min_colL,bwperim(mask_L1));
    if min_distL > 75 || min_rowL < 150
      left_mask = mask_L2 | mask_L1;
      left_mask(sub2ind(size(img),x,y)) = 1;
      mask_L2 = left_mask;
    else
      mask_L2(min_rowL:end,min_colL) = 1;
    end  

    [rowR2,colR2] = ind2sub(size(img),find(mask_R2));
    [min_rowR,min_indR] = max(rowR2);
    min_colR = colR2(min_indR);
    mask_R1(1:min_rowR,:) = 0;

    [x,y,min_distR] = minPointObjectDist(min_rowR,min_colR,bwperim(mask_R1));
    if min_distR > 75 || min_rowR < 150
      right_mask = mask_R2 | mask_R1;
      right_mask(sub2ind(size(img),x,y)) = 1;
      mask_R2 = right_mask;
    else
      mask_R2(min_rowR:end,min_colR) = 1;
    end

    mask_end = mask | mask_L2 | mask_R2;
  end
  
  mask_filled = imfill(mask_end,'holes');
  mask_last = mask_filled - mask_end;
  mask_last = bwareafilt(logical(mask_last),1,'largest',4);
  mask_last = bwconvhull(mask_last);

  %% Get Inner Well Image
  stats = regionprops(mask_last,'BoundingBox');
  left = ceil(stats.BoundingBox(1));
  top = ceil(stats.BoundingBox(2));
  width = min([500-left, stats.BoundingBox(3)]);
  height = min([500-top, stats.BoundingBox(4)]);

  cropped_img = img .* uint8(mask_last);
  cropped_img = cropped_img(top:top+height,left:left+width);
  
  imshow(cropped_img); title('Cropped Image'); pause(2)

  out_name = fullfile(diroutput,file_name);
  imwrite(mask_last, char(out_name));
end
fprintf('\nEnd of Program\n');