clc
clear
close all

%%% here we use the affine matrix to describe the conjugate relationship
%%% between the spatial light modulator and camera
savepath = [ 'rawdata\CCD2DMD\calibrated_CCD2DMD_conjugate_coordinates\'];
if ~exist(savepath)
    mkdir(savepath);
end
%%%------------read the captured calibration image -------------------
img1 = im2double(imread([ 'rawdata\CCD2DMD\captured_calibration_image.tif']));

img1 = mat2gray(img1);
[img_r,img_c] = size(img1);

%%% extract the corner points
img1_points = img1(:,1:img_c);
[CORNERS,~,~] = detectCheckerboardPoints(img1_points);
figure,imshow(img1,[]);
hold on;plot(CORNERS(:,1),CORNERS(:,2),'b.','MarkerSize',8);

%%% sort the corner ---------------------
corner_c = sortrows(CORNERS,1);
ncol = 15;
mrow = 15;
for ii = 1:ncol
    a1 = sortrows(corner_c((ii-1)*mrow+1:ii*mrow,:),2);
    corner_c((ii-1)*mrow+1:ii*mrow,:) =a1 ;
end
corner_row = corner_c(:,2);
corner_col = corner_c(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------------------------------
figure,imshow(img1,[]);
hold on;plot(corner_col,corner_row,'r*');

%%%%%%%%%%%%-----------read the projected calibration pattern---------------
pattern=im2double(imread('rawdata\CCD2DMD\projected_calibration_pattern.png'));

%%%%-----------generate the corner coordinates on the pattern projected by the DMD -----------------
tt = 120:60:960;
tt = sort(tt,'ascend');
k1 = 1;
points1 = zeros(225,2);
for nn = 512.5:64:1408.5 %% column coordinate
    %% row coordinate.  During the experiment, we rotated the camera by 180 degrees, which caused the captured calibration image to appear inverted relative to the displayed calibration image.
    for mm = 988.5:-64:92.5 
        points1(k1,1) = mm;
        points1(k1,2) = nn;
        k1 = k1+1;
    end
end

prow = points1(:,1);
pcol = points1(:,2);

figure,imshow(pattern);
hold on;plot(pcol(:),prow(:),'r.','MarkerSize',28);

%%%-------------construct matrix for determining the affine matrix --------------------
A = [pcol,prow];
B = [corner_col,corner_row];
N2 = [A,ones(size(A,1),1)];       
M2 = [B,ones(size(B,1),1)];        

 %%%%-----------calculate the affine matrix --------------------
aM =M2\N2;                

%%%---------------------calculate the reprojection error  ---------------
N2_new = M2*aM;
delta_N = (N2_new - N2);
%%% show the reprojection error-------------------
figure,plot(delta_N(:,1),delta_N(:,2),'r*');xlabel('x (Pixel)');ylabel('y (Pixel)');title('Reprojection error');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% size of the captured image
image_row=512;
image_column=512;

%%%% The images used for calibration are full - size images captured by the camera. 
%%%% However, during the actual measurement, the images captured by the camera are not full - size but cropped ones. 
%%%% Therefore, a coordinate transformation is required here.
[CCD_col,CCD_row] = meshgrid(img_c/2-image_column/2+1:img_c/2+image_column/2,img_r/2-image_row/2+1:img_r/2+image_row/2);

CCD_coordinate = [CCD_col(:),CCD_row(:)];
CCD_coordinate = [CCD_coordinate,ones(size(CCD_coordinate,1),1)];

%%% Based on the calibrated affine matrix, calculate the coordinates of the conjugate points of each pixel in the camera within the coordinate system of the spatial light modulator.
CCD_2_DMD = CCD_coordinate*aM;
figure,plot(CCD_2_DMD(:,1),CCD_2_DMD(:,2),'r*');

CCD2DMD_col = CCD_2_DMD(:,1) ;  %%% 列坐标
CCD2DMD_row = CCD_2_DMD(:,2) ;  %%% 行坐标

%%% size of the DMD  (unit: pixel)
DMD_ncol = 1920;
DMD_nrow = 1080;

%%% size of the projected pattern (unit: pixel)
pattern_nrow_DMD=1024;
pattern_ncol_DMD=1024;

%%%% The image load on the DMD are full size. 
%%%% However, the size of the pattern used for measurement are not full - size. 
%%%% Therefore, a coordinate transformation is required here.
CCD2DMD_sub_col0 = CCD2DMD_col - (DMD_ncol - pattern_ncol_DMD)/2;
CCD2DMD_sub_row0 = CCD2DMD_row - (DMD_nrow - pattern_nrow_DMD)/2;

%% To shorten the reconstruction time of the dual image, we reduced its size. 
%% Meanwhile, the size of the pattern required for single - pixel measurement was also correspondingly decreased. 
%% However, to avoid affecting the field-of-view measurement, we adopted the method of periodic extension to expand the pattern periodically.
%% Here, the periodic extension factor is set to 8, which means that the pattern will be extended 8 times in the horizontal and vertical directions respectively. 
%% Assuming that the size of the original pattern is 128×128 pixels, after the extension, the size of the pattern will become 1024×1024 pixels.
for block_num = [8]  %% periodic extension factor 
    
    %%% Calculate the size of the dual image based on the extension factor and the size of the actually projected pattern.
    block_size = pattern_nrow_DMD/block_num;
    
    %Transform the coordinates of the conjugate points of the camera on the spatial light modulator to the coordinate system of the dual image.
    CCD2DMD_sub_col = mod(CCD2DMD_sub_col0,block_size);
    CCD2DMD_sub_row = mod(CCD2DMD_sub_row0,block_size);
    
    %%% reshape and save the coordinates of the conjugate points
    CCD2DMD_sub_col_2D = reshape(CCD2DMD_sub_col,image_row,image_column);
    CCD2DMD_sub_row_2D = reshape(CCD2DMD_sub_row,image_row,image_column);
    
    D_row = zeros(image_row,image_column);
    D_col = zeros(image_row,image_column);
    for jj = 1:image_row
        for kk = 1:image_column
            Drow = CCD2DMD_sub_row_2D(jj,kk);
            Dcol = CCD2DMD_sub_col_2D(jj,kk);
            D_row(jj,kk) = Drow;
            D_col(jj,kk) = Dcol;
        end
    end
    
    save([savepath 'D_cor_s_0_block_Drow_',num2str(block_num),'.mat'],'D_row');
    save([savepath 'D_cor_s_0_block_Dcol_',num2str(block_num),'.mat'],'D_col');
end
