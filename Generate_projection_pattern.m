
clc
close all
clear

T_Pixel =1024;
%% size of the DMD
DMD_x=1920;
DMD_y=1080;
f_sp = 76;
l_sp = 99;
pattern_path = 'projection_pattern\';

if ~exist(pattern_path)
    mkdir(pattern_path);
end

Block = 8;
npixel = T_Pixel/Block;

mRow = T_Pixel/Block;
nCol  =T_Pixel/Block;
SamplingPath = 'circular';

%% Parameters
nStepPS = 3;
Phaseshift = 90;
Amplitude = 1;
SpectralCoverage = 1;
pattern_new1 = zeros(DMD_y,DMD_x);

[fxMat, fyMat] = meshgrid([0:1:nCol-1]/nCol, [0:1:mRow-1]/mRow);
fxMat = fftshift(fxMat);
fyMat = fftshift(fyMat);

OrderMat = getOrderMat(mRow, nCol, SamplingPath);
[nCoeft,tmp] = size(OrderMat);
nCoeft = round(nCoeft * SpectralCoverage);

InitPhaseArr = getInitPhaseArr(nStepPS, Phaseshift);
IntensityMat = zeros(mRow, nCol, nStepPS);

RealFourierCoeftList = getRealFourierCoeftList(mRow, nCol);



%% Main loop for generate time-varying patterns illumiantion
ii=1;
for iCoeft = f_sp:l_sp
    
    count= iCoeft;
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
    iRow = OrderMat(iCoeft,1);
    jCol = OrderMat(iCoeft,2);
    
    fx = fxMat(iRow,jCol);
    fy = fyMat(iRow,jCol);
    
    IsRealCoeft = existVectorInMat( [iRow jCol], RealFourierCoeftList );
    
    for iStep = 1:nStepPS
        if IsRealCoeft == 1 && iStep > 2
            if nStepPS == 3
                IntensityMat(iRow,jCol,iStep) = IntensityMat(iRow,jCol,2);
            end
            if nStepPS == 4
                IntensityMat(iRow,jCol,iStep) = 0;
            end
            continue;
        end
        
        [ sub_Pattern ] = getFourierPattern( Amplitude, mRow, nCol, fx, fy, InitPhaseArr(iStep) );    %
        
        %%Periodically extended pattern.
        extended_pattern = repmat(sub_Pattern,[8,8]);
        projected_pattern = zeros(DMD_y,DMD_x);
        projected_pattern(DMD_y/2-T_Pixel/2+1:DMD_y/2+T_Pixel/2,DMD_x/2-T_Pixel/2+1:DMD_x/2+T_Pixel/2) = extended_pattern;
        imwrite(projected_pattern,[pattern_path num2str(picnum) '.png']);
        picnum=picnum+1;
    end
    
end

