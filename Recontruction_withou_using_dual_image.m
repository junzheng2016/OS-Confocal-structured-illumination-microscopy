clc
clear
close all

T_Pixel = 1024; %%% The size of the projected structured pattern. (unit: pixel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% In the confocal structured illumination microscopy system, each pixel of the camera is regarded as an independent single - pixel detector.
%%%%%%%%%% The dual image is reconstructed using the Fourier single - pixel imaging method, with the core objective of effectively separating the conjugate signal from the non - conjugate signal noise.
%%%%%%%%%% After successfully separating the signals and noise, a digital pinhole is constructed based on the conjugate relationship between the spatial light modulator and the camera.
%%%%%%%%%% This digital pinhole is generated based on the physical characteristics and interrelationships of the components in the optical system, and it has clear optical significance and practical application value.
%%%%%%%%%% Finally, the generated digital pinhole is used to accurately extract the conjugate signal from the reconstructed dual image.
%%%%%%%%%% By rearranging the extracted signal according to the pixel coordinate of the camera pixel, we can reconstruct the confocal image.

%%%% The serial numbers of the spatial frequencies corresponding to the Fourier patterns used in the measurement.
f_sp = 76;   % the first one
l_sp = 99;   % the last one
Block =[8];    % the number of block

%%% Load the coordinates of the conjugate mapping points corresponding to each pixel of the camera in the pixel coordinate system of the spatial light modulator.
load  ([ 'rawdata\CCD2DMD\calibrated_CCD2DMD_conjugate_coordinates\D_cor_s_0_block_Drow_' num2str(Block) '.mat'],'D_row');
load  ([ 'rawdata\CCD2DMD\calibrated_CCD2DMD_conjugate_coordinates\D_cor_s_0_block_Dcol_' num2str(Block) '.mat'],'D_col');
irow_0 = D_row(:);   % row coordinate;
jcol_0 = D_col(:);   % column coordinate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% create the fold to save the result
result_savepath='rawdata\tz\pixel_1024_1024_block_8_8_coeft_76_99\result\OS_CSIM_wihout_using_dual_image\';

if ~exist(result_savepath)
    mkdir(result_savepath);
end

%%%%% result with normalization
result_savepath_nor=[ result_savepath 'normalized\'];
if ~exist(result_savepath_nor)
    mkdir(result_savepath_nor);
end

%%%%% result without normaliztion
result_savepath_unor=[ result_savepath 'unnormalized\'];
if ~exist(result_savepath_unor)
    mkdir(result_savepath_unor);
end

for depth = 0:100  %%% The depth range of the sample to be measured.
    tic
    depth
    
    nStepPS = 3;             % Phase shift steps
    Phaseshift = 120;   %% phase shift. (unit: degree)
    
    count= f_sp;
    picnum = 1;
    
    if nStepPS == 4
        if count==1
            picnum=1;
        else
            picnum= (count-1)*nStepPS-2+1;
        end
    end
    
    if nStepPS == 3
        if count == 1
            picnum = 1;
        else
            picnum = (count-1)*nStepPS;
        end
    end
    
    path = strcat('rawdata\tz\pixel_1024_1024_block_8_8_coeft_76_99\depth', num2str(depth,'%4.2f'), 'um_binary_16bit\');
    
    img0=im2double(imread([path 'Image_',num2str(picnum),'.png']));
    
    [picsizex,picsizey]=size(img0);         %   obtain image size
    
    npixel = T_Pixel/Block;
    
    mRow = npixel;
    nCol  =npixel;
    SamplingPath = 'circular';
    
    OrderMat = getOrderMat(mRow, nCol, SamplingPath);                              % generate sampling path in Fourier domain生成采样的路径的坐标

    
    InitPhaseArr = getInitPhaseArr(nStepPS, Phaseshift);                       %返回相位值
    RealFourierCoeftList = getRealFourierCoeftList(mRow, nCol);                  %找出实值傅里叶系数的频率所在坐标的位置
    
    [fxMat, fyMat] = meshgrid([0:1:nCol-1]/nCol, [0:1:mRow-1]/mRow);           % generate coordinates in Fourier domain (not neccessary) 生成每行每列的频率矩阵
    fxMat = fftshift(fxMat);
    fyMat = fftshift(fyMat);
    
    
    %%%%   Define a matrix to store the Fourier coefficients calculated by taking each pixel point of the camera as a single-pixel detector.
    AllCoeft = zeros(picsizex*picsizey,l_sp-f_sp+1);
    pic = zeros(picsizex,picsizey,nStepPS);
    
    ik = 1;
    interval = 1;
    for count= f_sp:interval:l_sp
        
        iRow = OrderMat(count,1);
        iCol = OrderMat(count,2);
        fxx(ik) = fxMat(iRow,iCol);
        fyy(ik) = fyMat(iRow,iCol);
        ik = ik + 1;
        
        if nStepPS == 4
            if count==1
                picnum=1;
            else
                picnum= (count-1)*nStepPS-2+1;
            end
        end
        
        if nStepPS == 3
            if count == 1
                picnum = 1;
            else
                picnum = (count-1)*nStepPS;
            end
        end
        
        IsRealCoeft = existVectorInMat( [iRow iCol], RealFourierCoeftList );  
        
        for iStep = 1:nStepPS  
            if IsRealCoeft == 1 && iStep > 2       
                if nStepPS == 4
                    pic(:,:,iStep)=zeros(picsizex,picsizey);
                end
                if nStepPS == 3
                    pic(:,:,iStep) = pic(:,:,2);
                end
                
                continue;
            end
            img1=(im2double(imread([path 'Image_' num2str(picnum,'%0.1d') '.png'])));
            
            picnum=picnum+1;
            
            pic(:,:,iStep)= (img1);%(position(1):position(2),position(3):position(4));
        end
        WF =  (pic(:,:,1) + pic(:,:,2) + pic(:,:,3))/3;
        
        if nStepPS == 4
            spec = (pic(:,:,1)-pic(:,:,2))+1i*(pic(:,:,3)-pic(:,:,4));
        end
        if nStepPS == 3
            spec = (2*pic(:,:,1)-pic(:,:,2)-pic(:,:,3))+sqrt(3)*1i*(pic(:,:,2)-pic(:,:,3));
        end
        
        AllCoeft (:,count-f_sp+1) = spec(:);
        
    end
    
    %%% calculate the argument of the Fourier coefficients.
    ang = angle(AllCoeft);
    %%% the modulus of the Fourier coefficients.
    modu = abs(AllCoeft);
    
    %%%  Define four matrices to store the coordinates of the four integer-pixel points (point1,point2,point3,point4) adjacent to the conjugate point.
    point1_row_id = zeros(picsizex*picsizey,1);point1_col_id = zeros(picsizex*picsizey,1);
    point2_row_id = zeros(picsizex*picsizey,1);point2_col_id = zeros(picsizex*picsizey,1);
    point3_row_id = zeros(picsizex*picsizey,1);point3_col_id = zeros(picsizex*picsizey,1);
    point4_row_id = zeros(picsizex*picsizey,1);point4_col_id = zeros(picsizex*picsizey,1);
    
    firow_id = zeros(picsizex*picsizey,1);fjcol_id = zeros(picsizex*picsizey,1);
    
    for num = 1:picsizex*picsizey
        
        irow = irow_0(num,1);  jcol = jcol_0(num,1);
        
        firow = fix(irow);
        fjcol = fix(jcol);
        
        firow_id(num,1) = firow;
        fjcol_id(num,1) = fjcol;
        
        %%% Determine the coordinates of the four integer-pixel points adjacent to the conjugate point.
        if firow > 0 && fjcol > 0 && firow < npixel && fjcol<npixel
            point1_row_id(num,1) = firow;       point1_col_id(num,1) = fjcol;
            point2_row_id(num,1) = firow+1;     point2_col_id(num,1) = fjcol;
            point3_row_id(num,1) = firow;       point3_col_id(num,1) = fjcol+1;
            point4_row_id(num,1) = firow+1;     point4_col_id(num,1) = fjcol+1;
        end
        
        if firow == 0 && fjcol == 0
            point1_row_id(num,1) = npixel;     point1_col_id(num,1) = npixel;
            point2_row_id(num,1) = 1;          point2_col_id(num,1) = npixel;
            point3_row_id(num,1) = npixel;     point3_col_id(num,1) = 1;
            point4_row_id(num,1) = 1;          point4_col_id(num,1) = 1;
        end
        
        if firow == 0 && fjcol ~= 0
            point1_row_id(num,1) = npixel;     point1_col_id(num,1) = fjcol;
            point2_row_id(num,1) = 1;          point2_col_id(num,1) = fjcol;
            point3_row_id(num,1) = npixel;     point3_col_id(num,1) = fjcol+1;
            point4_row_id(num,1) = 1;          point4_col_id(num,1) = fjcol+1;
        end
        
        if firow ~=0 && fjcol == 0
            point1_row_id(num,1) = firow;      point1_col_id(num,1) = npixel;
            point2_row_id(num,1) = firow+1;    point2_col_id(num,1) = npixel;
            point3_row_id(num,1) = firow;      point3_col_id(num,1) = 1;
            point4_row_id(num,1) = firow+1;    point4_col_id(num,1) = 1;
        end
        
    end
    
    %%% Subtract 1 from the coordinates because the counting starts from 0.
    point1_row_id = point1_row_id-1; point1_col_id = point1_col_id -1;
    point2_row_id = point2_row_id-1; point2_col_id = point2_col_id -1;
    point3_row_id = point3_row_id-1; point3_col_id = point3_col_id -1;
    point4_row_id = point4_row_id-1; point4_col_id = point4_col_id -1;
    
    %%%% Calculate the gray values of the four neighboring integer-pixel points.
    chazhi1 = sum(cos(2*pi*(repmat(fxx,picsizex*picsizey,1).*repmat(point1_col_id,1,l_sp-f_sp+1)+repmat(fyy,picsizex*picsizey,1).*repmat(point1_row_id,1,l_sp-f_sp+1))+ang).*modu,2);
    chazhi2 = sum(cos(2*pi*(repmat(fxx,picsizex*picsizey,1).*repmat(point2_col_id,1,l_sp-f_sp+1)+repmat(fyy,picsizex*picsizey,1).*repmat(point2_row_id,1,l_sp-f_sp+1))+ang).*modu,2);
    chazhi3 = sum(cos(2*pi*(repmat(fxx,picsizex*picsizey,1).*repmat(point3_col_id,1,l_sp-f_sp+1)+repmat(fyy,picsizex*picsizey,1).*repmat(point3_row_id,1,l_sp-f_sp+1))+ang).*modu,2);
    chazhi4 = sum(cos(2*pi*(repmat(fxx,picsizex*picsizey,1).*repmat(point4_col_id,1,l_sp-f_sp+1)+repmat(fyy,picsizex*picsizey,1).*repmat(point4_row_id,1,l_sp-f_sp+1))+ang).*modu,2);
    
    %%%% reshape the results
    chazhi01 = reshape(chazhi1,picsizex,picsizey);%
    chazhi02 = reshape(chazhi2,picsizex,picsizey);%
    chazhi03 = reshape(chazhi3,picsizex,picsizey);%
    chazhi04 = reshape(chazhi4,picsizex,picsizey);%
    
    
    %%% Calculate the gray value of the conjugate point using the bilinear interpolation method based on the gray values of the four neighboring integer-pixel points.
    pp = irow_0 - firow_id;
    qq = jcol_0 - fjcol_id;
    pp = reshape(pp,picsizex,picsizey);
    qq = reshape(qq,picsizex,picsizey);
    
    result_img = (1-pp).*(1-qq).*chazhi01 + pp.*(1-qq).*chazhi02 + (1-pp).*qq.*chazhi03 + pp.*qq.*chazhi04;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    result_img_unnormalized = result_img;
    result_img_neg_set_zeros = result_img_unnormalized;
    result_img_neg_set_zeros(result_img_neg_set_zeros<0)=0;   %%%%% Set the negative numbers in the result to zero.
    result_img_normalized = mat2gray(result_img_neg_set_zeros);  %% Normalize the result.
    toc
%     figure,imshow(result_img_normalized);
    %%%% save the results
    save([result_savepath_unor   'depth_', num2str(depth,'%4.2f'), 'um_unnormalize_sp_', num2str(f_sp),'_', num2str(l_sp),'_block_',num2str(Block),'x',num2str(Block),'.mat'],'result_img_unnormalized');

    imwrite(result_img_normalized,[result_savepath_nor  'depth_', num2str(depth,'%4.2f'), 'um_normalize_sp_',num2str(f_sp),'_',num2str(l_sp),'_block_',num2str(Block),'x',num2str(Block), '.bmp']);
    
end








