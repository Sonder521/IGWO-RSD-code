clear;
clc;
close all;

load('RSD.mat')
load('RSD2017.mat')
R = abs(RSD2017nian(:,1));
S = abs(RSD2017nian(:,2));
D = abs(RSD2017nian(:,3));

N = length(R);
E = zeros(N,1);
%求决策属性
for i = 1:N
    if R(i)<=1&&R(i)>=0.99
        E(i)=E(i)+0;
    elseif R(i)<0.99&&R(i)>=0.98
        E(i)=E(i)+1;
    elseif R(i)<0.98&&R(i)>=0.96
        E(i)=E(i)+2;
    elseif R(i)<0.96&&R(i)>=0.92
        E(i)=E(i)+3;
    elseif R(i)<0.92
        E(i)=E(i)+4;  
    end
    
    if S(i)<=0.1&&S(i)>=0
        E(i)=E(i)+0;
    elseif S(i)<=0.2&&S(i)>0.1
        E(i)=E(i)+1;
    elseif S(i)<=0.4&&S(i)>0.2
        E(i)=E(i)+2;
    elseif S(i)<=0.6&&S(i)>0.4
        E(i)=E(i)+3;
    elseif S(i)>0.6
        E(i)=E(i)+4;  
    end
    
    %暂不清楚D的决策向量应该如何赋值
    if D(i)<=0.05&&D(i)>=0
        E(i)=E(i)+0;
    elseif D(i)<=0.1&&D(i)>0.05
        E(i)=E(i)+1;
    elseif D(i)<=0.2&&D(i)>0.1
        E(i)=E(i)+2;
    elseif D(i)<=0.3&&D(i)>0.2
        E(i)=E(i)+3;
    elseif D(i)>0.3
        E(i)=E(i)+4;  
    end
    
end

%对条件属性进行归一化
% R0 = zeros(N,1);
% S0 = zeros(N,1);
% D0 = zeros(N,1);
% for i = 1:N
%     R0(i) = (R(i)-min(R))/(max(R)-min(R));
%     S0(i) = (S(i)-min(S))/(max(S)-min(S));
%     D0(i) = (D(i)-min(D))/(max(D)-min(D));  
% end

%RSDE写入一个矩阵

lidar_rsd = [R,S,D,E];
save('lidar.mat','lidar_rsd')
fprintf('写入完成，请返回主函数')


