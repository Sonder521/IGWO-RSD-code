% ���ڼ��㣨����ֲڼ���������������������Ҫ�ȵ�����ڳ���
% �����ӳ�����weigthD������Ȩ�أ���getPosSet����������reduceSet��Լ����㣩
%
%   ����˵���� �����һ�������жϣ��ж��������Ƿ����ĳ�о�����ͬ����
%             �����˺ܶ�����˵���� flag2 δ��������
%   �������ڣ�2015.04.06 �������ϴη������Ѿ�ʱ��һ���˰������˰���
%
% Made by suozi 20140428
% QQ��379786867

clc;
clear;
close all;

load('breast.mat') % ���ص�test.mat����Ϊ����ֲڼ�ʹ�õ���ֵ�;���ϵͳ����
load('lidar.mat')
%load('test.txt') % ͬѧ���ڼ������ݵ�ʱ��Ҳ����ֱ�Ӽ���txt�������ݣ�����excel��ʽ���ļ�
% ���һ��Ϊ��������
lammda=3; %����뾶������� delta=std��dataArray��/lammda
% ����lammdaȡֵ������0.5~1.5֮�䣬���̫���������������������̫С������򱨴�
% ��������ڰ������������Ƚ϶ࣨ��ʮ���ϣ��������lammda=2~4
sig_ctrl=0.001; %��Ҫ�����޵Ŀ��Ʋ�����ȡ�ӽ�0����

b = length(lidar_rsd);

c = randi(b,1,6000);
lidar_rsd1 = zeros(length(c),4);
for i = 1:length(c)
    lidar_rsd1(i,:) = lidar_rsd(c(i),:);
end
redSet = reduceSet(lidar_rsd1,lammda,sig_ctrl) %����Լ�򼯺�
weight = weightD(lidar_rsd1,lammda) %����Ȩ��
% redSet = reduceSet(lidar_rsd,lammda,sig_ctrl) %����Լ�򼯺�
% weight = weightD(lidar_rsd,lammda) %����Ȩ��


