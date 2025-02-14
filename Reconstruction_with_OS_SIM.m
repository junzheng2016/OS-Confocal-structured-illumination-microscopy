
clc
close all
clear

%%%% The serial numbers of the spatial frequencies corresponding to the Fourier patterns used in the measurement.
f_sp = 76;   % the first one
l_sp = 99;   % the last one

nStep=3;
phase_shift=2*pi/nStep;
pic_size = 512;

result_savepath=[ 'rawdata\tz\pixel_1024_1024_block_8_8_coeft_76_99\result\OS_SIM\result_sp_',num2str(f_sp),'_',num2str(l_sp),'\'];
if ~exist(result_savepath)
    mkdir(result_savepath);
end
path_nor = [result_savepath 'normalized\'];
if~exist(path_nor)
    mkdir(path_nor);
end

path_unor = [result_savepath 'unormalized\'];
if~exist(path_unor)
    mkdir(path_unor);
end

for depth =0:100
    
    path = strcat('rawdata\tz\pixel_1024_1024_block_8_8_coeft_76_99\depth', num2str(depth,'%4.2f'), 'um_binary_16bit\');

    if nStep == 4
        ps = [1,3,2,4];
    end
    if nStep == 3
        ps = [1,2,3];
    end
    
    result_avg = zeros(pic_size,pic_size);
    
    num_k = 0;
    for j = f_sp:l_sp
        
        %%% define three matries
        rsin=zeros(pic_size,pic_size);
        rcos=zeros(pic_size,pic_size);
        img_avg = zeros(pic_size,pic_size);
        
        for i=1:nStep
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read the image  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nStep == 4
                img_filename = strcat(path, 'Image_',num2str(ps(i)+(j-1)*nStep -2), '.png');
            end
            
            if nStep == 3
                img_filename = strcat(path, 'Image_',num2str(ps(i)+(j-1)*nStep-1), '.png');
            end
            
            img = im2double(imread(img_filename));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            rsin = img.*sin(phase_shift*(i-1))+rsin;
            rcos = img.*cos(phase_shift*(i-1))+rcos;
            img_avg = img_avg + img;
        end
        img_avg = img_avg/nStep;
        img_avg_nor = mat2gray(img_avg);
        
        %%%% calculate the optical sectioning result
        result = (rsin.^2+rcos.^2).^0.5;
        
        
        result_avg = result_avg + result;
        num_k = num_k + 1;
        
    end
    
    %%%%  average the optical sectioning result to reduce the noise
    result_avg = (result_avg/num_k);
    norresult_avg = mat2gray(result_avg);
    
    %%% save the result
    
    %%% unormalized result
    save([path_unor 'depth_',num2str(depth),'um_result_avg_' num2str(f_sp) '-' num2str(l_sp) '.mat'],'result_avg');
    %%%% normalized result
    imwrite(norresult_avg,[path_nor 'depth_',num2str(depth),'um_result_avg_' num2str(f_sp) '-' num2str(l_sp) '.bmp'])
    
end

