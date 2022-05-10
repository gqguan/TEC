%% 案例研究：绘制给定TEC在不同冷热侧温差下的性能
%
% by Dr. Guan Guoqiang @ SCUT on 2022-05-10

% 载入TEC参数
if exist('TEC','var') % 使用工作空间中的TEC参数
    opts = [1,0];
    fprintf('使用当前工作空间中的TEC参数\n')
else % 存盘文件中的TEC参数
    load('TEC_Params.mat','TEC_Params')
%     opts = [1,0]; TEC = TEC_Params.TEC(4);
%     opts = [1,0]; TEC = TEC_Params.TEC(6);
%     opts = [1,0]; TEC = TEC_Params.TEC(8);
%     opts = [1,0]; TEC = TEC_Params.TEC(13);
    opts = [1,0]; TEC = TEC_Params.TEC(15);
%     opts = [1,0]; TEC = TEC_Params.TEC(18);
end
% 初值
nGrid = 50;
opStr = 'max.All';
COP2 = zeros(1,nGrid);
Q2 = zeros(1,nGrid);
% 冷热侧温差范围
dT = linspace(1,50,nGrid);
Tc = 273.15+[30,20,10];
% 计算TEC特性
for j = 1:length(Tc)
    for i = 1:nGrid
        out = TE_ShowPerformance(Tc(j)+dT(i),Tc(j),TEC,opts,opStr);
        COP2(i) = out.MaxCOP2.Value;
        Q2(i) = out.MaxCOP2.Q2;
    end
    figure(1)
    plot(dT,Q2) % 考察不同T1对最高COP下的吸热量影响
    xlabel('T_h-T_c [K]')
    ylabel('Q_2@max.COP_2')
end
% 输出
figure(2)
yyaxis left
plot(dT,COP2)
xlabel('T_h-T_c [K]')
ylabel('max.COP_2')
yyaxis right
plot(dT,Q2)
ylabel('Q_2@max.COP_2 [W]')