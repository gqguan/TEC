%% 案例研究：绘制给定TEC在不同冷热侧温差下的性能
%
% by Dr. Guan Guoqiang @ SCUT on 2022-05-10

% 载入TEC参数
if exist('TEC','var') % 使用工作空间中的TEC参数
    opts = [0,0];
    fprintf('使用当前工作空间中的TEC参数\n')
else % 存盘文件中的TEC参数
    load('TEC_Params.mat','TEC_Params')
    iTEC = 15;
    TEC = TEC_Params.TEC(iTEC);
    opts = [TEC_Params.opt1(iTEC),TEC_Params.opt2(iTEC)];
%     opts = [1,0]; TEC = TEC_Params.TEC(4);
%     opts = [1,0]; TEC = TEC_Params.TEC(6);
%     opts = [1,0]; TEC = TEC_Params.TEC(8);
%     opts = [1,0]; TEC = TEC_Params.TEC(13);
%     opts = [0,1]; TEC = TEC_Params.TEC(14);
%     opts = [1,0]; TEC = TEC_Params.TEC(15);
%     opts = [1,0]; TEC = TEC_Params.TEC(18);
end
% 初值
nGrid = 50;
opStr = 'max.All';
COP2 = zeros(1,nGrid);
Q2 = zeros(1,nGrid);
maxQ2 = zeros(1,nGrid);
% 冷热侧温差范围
dT = linspace(1,50,nGrid);
Tc = 273.15+[30,20,10];
% 计算TEC特性
fig1 = figure('Units','pixels','Position',[100 100 1080 480]);
ax1 = axes;
ax2 = axes;
for j = 1:length(Tc)
    for i = 1:nGrid
        out = TE_ShowPerformance(Tc(j)+dT(i),Tc(j),TEC,opts,opStr);
        COP2(i) = out.MaxCOP2.Value;
        Q2(i) = out.MaxCOP2.Q2;
        maxQ2(i) = out.MaxQ2.Value;
    end
    subplot(1,2,1,ax1);
    p1 = plot(ax1,dT,maxQ2); % 考察不同T1对最高COP下的吸热量影响
    p1.DisplayName = sprintf('Tc = %.4g K',Tc(j));
    xlabel('T_h-T_c [K]')
    ylabel('max.Q_2')
    hold on
    subplot(1,2,2,ax2);
    yyaxis left
    p2left = plot(ax2,dT,COP2);
    p2left.DisplayName = sprintf('Tc = %.4g K',Tc(j)); 
    xlabel('T_h-T_c [K]')
    ylabel('max.COP_2')
    yyaxis right
    p2right = plot(ax2,dT,Q2);
    p2right.DisplayName = sprintf('Tc = %.4g K',Tc(j)); 
    ylabel('Q_2@max.COP_2 [W]')
    hold on
end
legend(ax1)
hold off
