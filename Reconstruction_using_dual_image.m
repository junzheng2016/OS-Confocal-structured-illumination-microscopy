clc
clear
close all

T_Pixel = 1024;%%% The size of the projected structured pattern. (unit: pixel)

%%%% The serial numbers of the spatial frequencies corresponding to the Fourier patterns used in the measurement.
f_sp = 76;   % the first one
l_sp = 99;   % the last one
Block =[8];    % the number of block

% %%% Load the coordinates of the conjugate mapping points corresponding to each pixel of the camera in the pixel coordinate system of the spatial light modulator.
% load  ([ 'rawdata\CCD2DMD\CCD2DMD_cord1\D_cor_s_0_block_Drow_' num2str(Block) '.mat'],'D_row');
% load  ([ 'rawdata\CCD2DMD\CCD2DMD_cord1\D_cor_s_0_block_Dcol_' num2str(Block) '.mat'],'D_col');
load  ([ 'rawdata\CCD2DMD\calibrated_CCD2DMD_conjugate_coordinates\D_cor_s_0_block_Drow_' num2str(Block) '.mat'],'D_row');
load  ([ 'rawdata\CCD2DMD\calibrated_CCD2DMD_conjugate_coordinates\D_cor_s_0_block_Dcol_' num2str(Block) '.mat'],'D_col');
irow_0 = D_row(:);   % row coordinate;
jcol_0 = D_col(:);   % column coordinate;


%%%%% create the fold to save the result
result_savepath='rawdata\tz\pixel_1024_1024_block_8_8_coeft_76_99\result\OS_CSIM_using_dual_image\';

if ~exist(result_savepath)
    mkdir(result_savepath);
end

%%%%% fold used to save result with normalization
result_savepath_nor=[ result_savepath 'normalized\'];
if ~exist(result_savepath_nor)
    mkdir(result_savepath_nor);
end

%%%%% fold used to save result without normaliztion
result_savepath_unor=[ result_savepath 'unnormalized\'];
if ~exist(result_savepath_unor)
    mkdir(result_savepath_unor);
end


for depth = 1:100
    depth
    
    path = strcat('rawdata\tz\pixel_1024_1024_block_8_8_coeft_76_99\depth', num2str(depth,'%4.2f'), 'um_binary_16bit\'); 
    img0=im2double(imread([path 'Image_225.png']));
    [picsizex,picsizey]=size(img0);         %   obtain image size

    npixel = T_Pixel/Block;
  
    mRow = npixel;
    nCol  =npixel;
    SamplingPath = 'circular';                                                 %%% sampling path
    nStepPS = 3;                                                               % number of phase shifting
    Phaseshift = 120;                                                           % phase shift
    
    OrderMat = getOrderMat(mRow, nCol, SamplingPath);                          % Generate the coordinates of the sampling path.
    
    InitPhaseArr = getInitPhaseArr(nStepPS, Phaseshift);                       % return the phase shift
    RealFourierCoeftList = getRealFourierCoeftList(mRow, nCol);                % Find the positions of the coordinates where the frequencies of the real-valued Fourier coefficients are located.
    
    result_img = zeros(picsizex,picsizey);                               %% matrix used to to store confocal image
    AllCoeft = zeros(picsizex*picsizey,l_sp-f_sp+1);                    %% matrix used to to store Fourier coefficients
    
    picnum = 1;
    pic = zeros(picsizex,picsizey,nStepPS);
    
    for count= f_sp:1:l_sp
        
        %%% Determine the numbers of the images to be read based on the numbers of the spatial frequencies.
        if nStepPS == 4
            if count==1
                picnum=1;
            else
                picnum= (count-1)*nStepPS-1;
            end
        end
        
        if nStepPS == 3
            if count == 1
                picnum = 1;
            else
                picnum = (count-1)*nStepPS;
            end
        end
        
        iRow = OrderMat(count,1);
        iCol = OrderMat(count,2);
        
        IsRealCoeft = existVectorInMat( [iRow iCol], RealFourierCoeftList );   
        
        %% read image
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
            img1=im2double(imread([path 'Image_' num2str(picnum,'%0.1d') '.png']));

            picnum=picnum+1;
            
            pic(:,:,iStep)= img1;
        end
        
        %%% use the four-step phase-shifting algorithm to calculate the Fourier
        %%% coefficient
        if nStepPS == 4
            fourier_coef = (pic(:,:,1)-pic(:,:,2))+1i*(pic(:,:,3)-pic(:,:,4));
        end
          %%% use the three-step phase-shifting algorithm to calculate the Fourier
        %%% coefficient
        if nStepPS == 3
            fourier_coef = (2*pic(:,:,1)-pic(:,:,2)-pic(:,:,3))+sqrt(3)*1i*(pic(:,:,2)-pic(:,:,3));
            fourier_coef = 2*fourier_coef/3;
        end
        
        AllCoeft (:,count-f_sp+1) = fourier_coef(:);
        
    end
    
    SpecMat = zeros(npixel,npixel);
    %%% Generate the index numbers of the measured Fourier coefficients in the frequency spectrum.
    Indx = sub2ind(size(SpecMat),OrderMat(f_sp:l_sp,1),OrderMat(f_sp:l_sp,2));
    
    %%% Reconstruct the dual images pixel by pixel and extract the
    %%% conjugate signals to assemble the confocal image
    for num = 1:picsizex*picsizey
        num
        SpecMat = zeros(npixel,npixel);
        SpecMat(Indx) = AllCoeft(num,:);   %%% assemble spectrum of the dual image

        fullSpec = completeSpec(SpecMat);   %%% Obtain the complete frequency spectrum using conjugate symmetry.
        dual_img =   real(ifft2(ifftshift(fullSpec)));   %%% reconstruct dual image by perform inverse Fourier transform
%         figure,imshow(dual_img,[]);
 
        %%%% extract conjugate signal from the dual image using bilinear
        %%%% interpolation algorithm
         irow = irow_0(num,1);  jcol = jcol_0(num,1);
        
        if jcol<=0
            jcol = jcol+npixel;
        end
        
        if jcol>=npixel
            jcol = jcol-npixel;
        end
        
        if irow<=0
            irow = irow+npixel;
        end
        
        if irow>=npixel
            irow = irow-npixel;
        end
        
        firow = fix(irow);
        fjcol = fix(jcol);
        %%% Determine the coordinates of the four integer pixel points adjacent to the conjugate point and use the bilinear interpolation method to calculate the conjugate signal.
        if firow > 0 && fjcol > 0 && firow < npixel && fjcol<npixel
            p = irow - firow;
            q = jcol - fjcol;
            confocal_value = (1-p)*(1-q)*dual_img(firow,fjcol)+p*(1-q)*dual_img(firow+1,fjcol)+(1-p)*q*dual_img(firow,fjcol+1)+p*q*dual_img(firow+1,fjcol+1);
        end
        
        if firow == 0 && fjcol == 0
            p = irow;
            q = jcol;
            confocal_value = (1-p)*(1-q)*dual_img(end,end)+p*(1-q)*dual_img(1,end)+(1-p)*q*dual_img(end,1)+p*q*dual_img(1,1);
        end
        
        if firow == 0 && fjcol ~= 0
            p = irow;
            q = jcol - fjcol;
            confocal_value = (1-p)*(1-q)*dual_img(end,fjcol)+p*(1-q)*dual_img(1,fjcol)+(1-p)*q*dual_img(end,fjcol+1)+p*q*dual_img(1,fjcol+1);
        end
        
        if firow ~=0 && fjcol == 0
            p = irow - firow;
            q = jcol;
            confocal_value = (1-p)*(1-q)*dual_img(firow,end)+p*(1-q)*dual_img(firow+1,end)+(1-p)*q*dual_img(firow,1)+p*q*dual_img(firow+1,1);
        end
        
        %%%%Place the extracted conjugate signals back to the coordinates where the camera pixels are located.
         [yy,xx] = ind2sub([picsizex,picsizey], num);
        result_img(yy,xx)= confocal_value;
  
    end
    
    result_img_unnormalized = result_img;
    result_img_neg_set_zeros = result_img_unnormalized;
    result_img_neg_set_zeros(result_img_neg_set_zeros<0)=0;   %%%%% Set the negative numbers in the result to zero.
    
    result_img_normalized = mat2gray(result_img_neg_set_zeros);  %% Normalize the result.
    toc
    
    %%%% save the results
    save([result_savepath_unor   'depth_', num2str(depth,'%4.2f'), 'um_unnormalize_sp_', num2str(f_sp),'_', num2str(l_sp),'_block_',num2str(Block),'x',num2str(Block),'.mat'],'result_img_unnormalized');
    imwrite(result_img_normalized,[result_savepath_nor  'depth_', num2str(depth,'%4.2f'), 'um_normalize_sp_',num2str(f_sp),'_',num2str(l_sp),'_block_',num2str(Block),'x',num2str(Block), '.bmp']);
 
end






