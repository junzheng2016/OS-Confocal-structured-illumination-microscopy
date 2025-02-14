%% GENERATE THE PATH MATRIX -> orderMat
% Rev 1: 2016/06/23 by Charles Cheung
% nPixel can be odd or even.
function OrderMat = getOrderMat(mRow, nCol, PathStr)
if mRow > nCol
    nPixel  = mRow;
else
    nPixel = nCol;
end

switch PathStr
    case 'spiral'
        PathMat = rot90(spiral(nPixel), 2);
    case 'circular'
        PathMat = CircularPath( nPixel, nPixel, floor(nPixel/2)+1, floor(nPixel/2)+1, ceil(ceil(nPixel/2)*sqrt(2)) );
    case 'diamond'
        PathMat = ZigzagFourier(nPixel);
end

nfc_s = 4;


[specx,specy] = meshgrid(1:mRow,1:nCol);
specx = specx-fix(mRow/2)-1;
specy = specy-fix(nCol/2)-1;
[theta,rho] = cart2pol(specx,specy);
%                             mask = 1;
mask =  rho>nfc_s;

temp= mask.*PathMat;

LeftMargin =  floor((nPixel - nCol)/2);                                %向下取整  leftmargin=左边距.floor?
RightMargin = floor((nPixel - nCol)/2);
TopMargin =  floor((nPixel - mRow)/2);
BottomMargin = floor((nPixel - mRow)/2);

PathMat = PathMat(TopMargin + 1 : TopMargin + mRow, ...
    LeftMargin + 1  : LeftMargin + nCol);

OrderMat = PathMat2OrderMat(PathMat, mRow, nCol, PathStr);
