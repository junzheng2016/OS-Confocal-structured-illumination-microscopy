function [ existence ] = existVectorInMat( vector, mat )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
mRow = size(mat, 1);                                     %提取实值（傅里叶系数点）矩阵的行数大小给mRow；size(A,1)表示返回A矩阵的行数
pz  = ones(mRow, 1) * vector - mat;                       %产生一个mRow一列的全1矩阵
existence = logical(size(find(sum(abs(pz)')' == 0), 1)); %logical将数值变为逻辑值 非零值为1，零值变0
end

