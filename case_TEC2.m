%% 案例研究：绘制给定TEC在不同冷热侧温度下的吸放热量及制冷系数曲面
%
% by Dr. Guan Guoqiang @ SCUT on 2022-04-20

%% 使用工作空间中的TEC参数
if exist('TEC','var') && exist('opts','var')
    fprintf('使用当前工作空间中的TEC参数\n')
else
    error('未在当前工作空间中找到变量TEC！')
end
%
nGrid = 50;
% 冷热侧温度范围
TCs = 273.15+linspace(15,35,nGrid);
THs = 273.15+linspace(40,60,nGrid);
% TEC输入电压
Vset = 9;
% 计算TEC吸放热量
Q1 = zeros(length(TCs),length(THs));
Q2 = zeros(size(Q1));
E = zeros(size(Q1));
COP2 = zeros(size(Q1));
for i = 1:length(TCs)
    for j = 1:length(THs)
        Tc = TCs(i);
        Th = THs(j);
        Q = TE_Heat(Th,Tc,TEC,opts(1),opts(2));
        Q1(i,j) = Q(1);
        Q2(i,j) = Q(2);
        E(i,j) = Q(1)-Q(2);
        COP2(i,j) = Q2(i,j)/E(i,j);
    end
end

%% 输出
[X,Y] = meshgrid(THs,TCs);
figure(1)
subplot(2,1,1)
surf(X,Y,COP2)
xlabel('T_H [K]')
ylabel('T_C [K]')
zlabel('COP_2')
subplot(2,1,2)
surf(X,Y,Q2)
xlabel('T_H [K]')
ylabel('T_C [K]')
zlabel('Q_2')
title(sprintf('U=%.1f V',Vset))