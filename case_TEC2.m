%% 案例研究：绘制给定TEC在不同冷热侧温度下的吸放热量曲面
%
% by Dr. Guan Guoqiang @ SCUT on 2022-04-20

%% 使用工作空间中的TEC参数
if exist('TEC','var')
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
for i = 1:length(TCs)
    for j = 1:length(THs)
        Tc = TCs(i);
        Th = THs(j);
        Q = TE_Heat(Th,Tc,TEC,0,1);
        Q1(i,j) = Q(1);
        Q2(i,j) = Q(2);
    end
end

%% 输出
[X,Y] = meshgrid(TCs,THs);
surf(X,Y,Q1)
xlabel('T_C [K]')
ylabel('T_H [K]')
zlabel('Q_1 [W]')
title(sprintf('U=%.1f V',Vset))