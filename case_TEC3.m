%% 案例研究：评估TEC_Params中的TEC参数，列出在不同冷热侧温差下的性能
%
% by Dr. Guan Guoqiang @ SCUT on 2022-05-11

%% 初始化
clear
nGrid = 50;
opStr = 'max.All';
COP2 = zeros(1,nGrid);
Q2 = zeros(1,nGrid);
maxQ2 = zeros(1,nGrid);
outTab = table;
% 冷热侧温差范围
dT = linspace(15,45,nGrid);
Tc = 273.15+[30,20,10]';
% 载入TEC参数
load('TEC_Params.mat','TEC_Params')

%% 计算各TEC特性
MaxCOP2Q2 = zeros(size(Tc));
for i = 1:height(TEC_Params)
    if ~isequal(TEC_Params.Status{i},'normal')
        continue
    end
    TEC = TEC_Params.TEC(i);
    opts = [TEC_Params.opt1(i),TEC_Params.opt2(i)];
    for j = 1:length(Tc)
        out(j,:) = arrayfun(@(x)TE_ShowPerformance(Tc(j)+x,Tc(j),TEC,opts,opStr),dT);
    end
    [MaxCOP2,iMaxCOP2] = max(arrayfun(@(x)x.MaxCOP2.Value,out),[],2);
    rawMaxCOP2Q2 = arrayfun(@(x)x.MaxCOP2.Q2,out);
    for k = 1:length(Tc)
        MaxCOP2Q2(k) = rawMaxCOP2Q2(1,iMaxCOP2(k));
    end
    THmaxCOP2 = Tc+dT(iMaxCOP2)';
    [MinCOP2,iMinCOP2] = min(arrayfun(@(x)x.MaxCOP2.Value,out),[],2);
    THminCOP2 = Tc+dT(iMinCOP2)';
    [MaxQ2,iMaxQ2] = max(arrayfun(@(x)x.MaxQ2.Value,out),[],2);
    THmaxQ2 = Tc+dT(iMaxQ2)';
    [MinQ2,iMinQ2] = min(arrayfun(@(x)x.MaxQ2.Value,out),[],2);
    THminQ2 = Tc+dT(iMinQ2)';
    ITEC(1:length(Tc),1) = i;
    outTab = [outTab;table(ITEC,Tc,MaxCOP2,THmaxCOP2,MaxCOP2Q2,MinCOP2,THminCOP2,MaxQ2,THmaxQ2,MinQ2,THminQ2)];
end
[MaxMinQ2,iMax] = max(outTab.MaxQ2);
fprintf('在TEC_Params中ITEC=%d的TEC参数在Tc=%.4g[K]可获得最大的最低吸热量Q2=%.4g[W]\n',...
    outTab.ITEC(iMax),outTab.Tc(iMax),MaxMinQ2)
disp(outTab(iMax,:))
opts = [TEC_Params.opt1(outTab.ITEC(iMax)),TEC_Params.opt2(outTab.ITEC(iMax))];
PlotResults(Tc,dT,TEC_Params.TEC(outTab.ITEC(iMax)),opts,opStr)
% for j = 1:length(Tc)
%     for i = 1:nGrid
%         out = TE_ShowPerformance(Tc(j)+dT(i),Tc(j),TEC,opts,opStr);
%         COP2(i) = out.MaxCOP2.Value;
%         Q2(i) = out.MaxCOP2.Q2;
%         maxQ2(i) = out.MaxQ2.Value;
%     end
% 
% end

%% 画图
function PlotResults(Tc,dT,TEC,opts,opStr)
    figure('Units','pixels','Position',[100 100 1080 480]);
    ax1 = axes;
    ax2 = axes;
    for j = 1:length(Tc)
        out = arrayfun(@(x)TE_ShowPerformance(Tc(j)+x,Tc(j),TEC,opts,opStr),dT);
        maxQ2 = arrayfun(@(x)x.MaxQ2.Value,out);
        maxCOP2 = arrayfun(@(x)x.MaxCOP2.Value,out);
        Q2 = arrayfun(@(x)x.MaxCOP2.Q2,out);
        subplot(1,2,1,ax1);
        p1 = plot(ax1,dT,maxQ2); % 考察不同T1对最高COP下的吸热量影响
        p1.DisplayName = sprintf('Tc = %.4g K',Tc);
        xlabel('T_h-T_c [K]')
        ylabel('max.Q_2')
        hold on
        subplot(1,2,2,ax2);
        yyaxis left
        p2left = plot(ax2,dT,maxCOP2);
        p2left.DisplayName = sprintf('Tc = %.4g K',Tc); 
        xlabel('T_h-T_c [K]')
        ylabel('max.COP_2')
        yyaxis right
        p2right = plot(ax2,dT,Q2);
        p2right.DisplayName = sprintf('Tc = %.4g K',Tc); 
        ylabel('Q_2@max.COP_2 [W]')
        hold on
        legend(ax1)
    end
    hold off
end
