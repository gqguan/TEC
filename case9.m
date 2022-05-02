%% 模拟一维集成热泵DCMD膜组件中的传热和传质现象
%
% by Dr. Guan Guoqiang @ SCUT on 2022/05/02

%% 初始化
clear
% 调用公用变量定义，其中包括DuctGeom（流道几何尺寸）、Stream（物料定义）、MembrProps（膜材料性质）
CommonDef
% 设定集成TEC多级SFMD系统的级数
NumStage = 1;
% 设定膜组件的热、冷侧进口温度和流率
T1 = 323.15; T2 = 288.15; % [K]
W1 = 1.5e-4; W2 = 1.5e-4; % [kg/s]
% 设定环境温度
TEXs = [298.15,298.15]; % [K]
% 膜组件热侧进料初始化
s1 = Stream;
s1.Temp = T1;
s1.MassFlow = W1;
% calculate the rest properties of feed-side influent
sIn(1) = DCMD_PackStream(s1);
% 膜组件冷侧进料初始化
s2 = Stream;
s2.Temp = T2;
s2.MassFlow = W2;
% calculate the rest properties of permeate-side influent
sIn(2) = DCMD_PackStream(s2);
% 设定各级膜材料特性
membrane = MembrProps;
% set properties for all TECs
load('TEC_Params.mat') % 载入已有的TEC计算参数
TECs(1:(NumStage+1)) = TEC_Params.TEC(3,1); % 注意按opt=0计算TEC的吸放热量

%% 计算集成热泵DCMD膜组件中的温度分布
[profile1,sOut1] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"countercurrent");
[profile2,sOut2] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"cocurrent");

%% 输出
DispResults(profile1,DuctGeom,membrane)
DispResults(profile2,DuctGeom,membrane)

function DispResults(profile,DuctGeom,membrane)
    persistent iLine
    if isempty(iLine)
        iLine = 0;
    end
    iLine = iLine+1;
    % 产水量
    WP = sum([profile.SM.MassFlow]); % [kg/s]
    JM = WP/membrane.Area*3600; % [kg/m2-h]
    fprintf('产水量为%.4e[kg/s]（渗透通量为%.3f[kg/m2-h]）\n', WP, JM)
    % TEC特性
    E = sum(cell2mat(cellfun(@(x)x(:,1)-x(:,2),profile.QTEC,'UniformOutput',false)),2);
    Q1 = sum(cell2mat(cellfun(@(x)x(:,1),profile.QTEC,'UniformOutput',false)),2);
    Q2 = sum(cell2mat(cellfun(@(x)x(:,2),profile.QTEC,'UniformOutput',false)),2);
    fprintf('DCMD膜组件料液侧加热TEC电功率为%.4f[W]：放热量%.4f[W]，从环境吸热量%.4f[W]\n', ...
        E(1), Q1(1), Q2(1))
    fprintf('DCMD膜组件渗透侧冷却TEC电功率为%.4f[W]：吸热量%.4f[W]，向环境放热量%.4f[W]\n', ...
        E(2), Q2(2), Q1(2))
    % 绘制温度侧型
    lineColor = {'#7E2F8E';'#A2142F';'#D95319';'#77AC30';'#4DBEEE';'#0072BD'};
    lineStyle = {'-';'--';':';'-.'};
    x = linspace(0,DuctGeom.Length,length(profile.QM));
    p = plot([x;x;x;x;x;x]',profile.T',lineStyle{iLine});
    for iPlot = 1:length(p)
        p(iPlot).Color = lineColor{iPlot};
    end
    if iLine <= 4
        hold on
    else
        hold off
        iLine = 0;
    end
end

% %% 计算SFMD性能
% % 定义优化目标函数为单位能耗
% opfun = @(x)interface_SFMD(x, NumStage, SInFeeds, SInPerms, Membranes, TECs);
% % 初值
% x0 = [W1 T1 I1 W2 T2 I2];
% % 变量边界
% lb = [3.5e-6, 40+273.15, 0, 3.5e-6, 5+273.15, 0];
% ub = [2.1e-3, 80+273.15, 8, 2.1e-3, 30+273.15, 8];
% % 求优
% opts = optimoptions('fmincon','Display','iter');
% x = fmincon(opfun,x0,[],[],[],[],lb,ub,[],opts);
% 
% %% 输出结果
% [SEC,out1,out2] = opfun(x);
% format short g
% fprintf('Specific energy consumption of %d-stage DCMD is %.4e J/kg (%.4g kWh/kg).\n', NumStage, SEC, SEC*2.778e-7);
% disp(out1)