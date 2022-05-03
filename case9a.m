%% 模拟一维集成热泵DCMD膜组件中的传热和传质现象：传统DCMD系统的单位能耗
% 当膜组件集成半导体热泵TEHP时，可载入相应的TEC参数，例如用TEC_Params.mat中的H28、H05等
% 当未集成TEHP时，可载入TEC_Params.mat中的H00（近似绝热的边界条件：无电功输入且导热系数很小）
%
% by Dr. Guan Guoqiang @ SCUT on 2022/05/02
%
% 在稳态操作条件下，料液加热所需热量Q1为使膜组件料液侧出口温度升高到膜组件进口温度；
% 渗透液冷却所需冷量Q2为膜组件渗透侧出口温度降低到膜组件进口温度。
% 由此，料液加热器能耗按电功100%转化为热量计算，即为E1 = Q1；
% 渗透液冷却器能耗按TEC制冷耗电量计算，即为E2 = Q2/COP2
%
% 对于最低单位能耗的全回流稳态操作，料液加热量Q1 = QM+WF*CP*(TMH-T0)，详细推导见笔记2022/5/3；
% 渗透液冷却吸热量Q2 = QM+WP*CP*(TMC-TP2)
%
% 为确定DCMD系统吸放热量Q1和Q2，需要获得跨膜传热量和渗透量（QM和WP）的基础上进一步确定膜面温度（TMH和TMC）
% 方法1：分别列出WP和QM的计算式联立求解得TMH和TMC；
% 方法2：对模拟结果中的温度分布（即T(3)和T(4)）求平均值

%% 初始化
clear
% 调用公用变量定义，其中包括DuctGeom（流道几何尺寸）、Stream（物料定义）、MembrProps（膜材料性质）
CommonDef
% 设定集成TEC多级SFMD系统的级数
NumStage = 1;
% 设定膜组件的热、冷侧进口温度和流率
T1 = 323.15; T2 = 288.15; % [K]
W1 = 1.5e-4; W2 = 1.5e-5; % [kg/s]
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
% 注意按opt1=0,opt2=1计算TEC的吸放热量
% opts = [0,1]; TECs(1:(NumStage+1)) = TEC_Params.TEC(18,1);
% opts = [1,0]; TECs(1:(NumStage+1)) = TEC_Params.TEC(14,1);
% opts = [0,0]; TECs(1:(NumStage+1)) = TEC_Params.TEC(4,1);
opts = [0,0]; TECs(1:(NumStage+1)) = TEC_Params.TEC(1,1); % 相当于未集成半导体热泵的DCMD膜组件


%% 计算集成热泵DCMD膜组件中的温度分布
[profile1,sOut1] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"countercurrent",opts);
% [profile2,sOut2] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"cocurrent",opts);
opStr = 'cooling';

%% DCMD系统单位能耗
% 计算稳态操作时料液放热量Q(1)和渗透液吸热量Q(2)
[Q,WP,~,TP1,TP2] = CalcHeat(profile1);
% 加热器功耗
E(1) = Q(1);
% 计算膜组件外置半导体制冷功耗
% opts = [0,1]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(18,1);
opts = [1,0]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(14,1);
% opts = [0,0]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(4,1);
[E(2),QTEC,QDiff] = CalcTECPower(opStr,Q(2),TEXs(1),mean([TP1,TP2]),exTECs(1),opts);
% 计算系统总能耗
SEC = sum(E)/WP/3600/1000; % [kWh/kg]

%% 输出
fprintf('DCMD系统在膜组件进口温度分别为%.2f[K]和%.2f[K]、流率分别为%.4g[kg/s]和%.4g[kg/s]\n',T1,T2,W1,W2)
fprintf('DCMD系统稳定操作需要加热量%.4g[W]和制冷量%.4g[W]\n',Q(1),Q(2))
fprintf('系统制冷采用TEC（编号%s）%s，计算热量为%.4g[W]，其偏差为%.4g[W]\n',TEC_Params.pid{14},opStr,QTEC,QDiff)
fprintf('DCMD系统产水率为%.4g[kg/h]，加热消耗能量%.4g[W]，冷却消耗能量%.4g[W]，单位能耗为%.4g[kWh/kg]\n',WP*3600,E(1),E(2),SEC)
figObj1 = DispResults(profile1,DuctGeom,membrane);
% figObj2 = DispResults(profile2,DuctGeom,membrane);

function p = DispResults(profile,DuctGeom,membrane)
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

function [Q,WP,QM,TP1,TP2] = CalcHeat(profile)
    T0 = 298.15; % 进料温度为环境温度[K]
    QM = sum(profile.QM); % 跨膜传热量[W]
    WP = sum(arrayfun(@(x)x.MassFlow,profile.SM)); % 跨膜渗透量[kg/s]
    TM = mean(profile.T(3:4,:),2);
    iStart = strfind(profile.Remarks,'：');
    switch profile.Remarks(iStart+1:end)
        case('cocurrent')
            TP1 = profile.S2(1).Temp;
            TP2 = profile.S2(end).Temp;
        case('countercurrent')
            TP1 = profile.S2(end).Temp;
            TP2 = profile.S2(1).Temp;
        otherwise
            error('CalcHeat()输入参数profile字段Remarks中无有效的流型信息')
    end
    cp1 = mean([profile.S1.SpecHeat]);
    cp2 = mean([profile.S2.SpecHeat]);
    Q(1) = QM+WP*cp1*(TM(1)-T0);
    Q(2) = QM+WP*cp2*(TM(2)-TP2);
end

function [E,Q,fval] = CalcTECPower(opStr,Q,Th,Tc,TEC,opts)
    switch opStr
        case('cooling')
            x0 = 2; lb = 0.1; ub = 12.5;
            solOpts = optimoptions(@lsqnonlin, 'Display', 'none');
            fun = @(x)GetTECHeat(x,opStr,Th,Tc,TEC,opts)-Q;
            x1 = lsqnonlin(fun,x0,lb,ub,solOpts);
            fval = fun(x1);
            [Q,E] = GetTECHeat(x1,opStr,Th,Tc,TEC,opts);
        case('heating')
        otherwise
    end
end

function [Qout,E] = GetTECHeat(var,opStr,Th,Tc,TEC,opts)
    switch opts(2)
        case(0)
            TEC.Current = var;
        case(1)
            TEC.Voltage = var;
        otherwise
            error('CalcTECHeat()输入参数opts有误！')
    end
    Q = TE_Heat(Th,Tc,TEC,opts(1),opts(2));
    E = Q(1)-Q(2);
    switch opStr
        case('cooling')
            Qout = Q(2);
        case('heating')
            Qout = Q(1);
        otherwise
            error('CalcTECHeat()输入参数opStr有误！')
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