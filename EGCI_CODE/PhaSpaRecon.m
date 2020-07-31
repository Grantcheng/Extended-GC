function [xn,L] = PhaSpaRecon(s,tau,m)
% 混沌序列的相空间重构 (phase space reconstruction)
% [xn,dn] = PhaSpaRecon(s,tau,m)
% 输入参数：    s          混沌序列
%               tau        重构时延
%               m          重构维数
% 输出参数：    xn         相空间中的点序列(每一列为相空间中一个点)
%               dn         一步预测的目标

len = length(s);
if (len-(m-1)*tau < 1)
    disp('err: delay time or the embedding dimension is too large!')
    xn = [];
    dn = [];
else
    xn = zeros(m,len-(m-1)*tau);
    L=len-(m-1)*tau;
    for i = 1:m
        %xn(i,:) = s(1+(i-1)*tau : len-(m-i)*tau);   % 相空间重构，每一列为一个点 
        xn(i,:) = s(len-(i-1)*tau :-1:1+(m-i)*tau);   % 相空间重构，每一列为一个点 
    end
    %xn=fliplr(xn);
end
