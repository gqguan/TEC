%% 多级集成热泵平板膜蒸馏单元的操作优化：一级SFMD最低单位能耗
%
% by Dr. Guan Guoqiang @ SCUT on 2020/11/25

%% 初始化
clear
% 调用公用变量定义，其中包括DuctGeom（流道几何尺寸）、Stream（物料定义）、MembrProps（膜材料性质）
CommonDef
% 设定集成TEC多级SFMD系统的级数
NumStage = 1;
% 设定膜组件的热、冷侧进口温度和流率
T1 = 323.15; T2 = 303.15;
W1 = 1.5e-4; W2 = 1.5e-4;
% 膜组件热侧进料初始化
SInFeed = Stream;
SInFeed.Temp = T1;
SInFeed.MassFlow = W1;
% calculate the rest properties of feed-side influent
SInFeed = DCMD_PackStream(SInFeed);
% 设定各级膜组件热侧进料
SInFeeds(1:NumStage) = SInFeed;
% 膜组件冷侧进料初始化
SInPerm = Stream;
SInPerm.Temp = T2;
SInPerm.MassFlow = W2;
% calculate the rest properties of permeate-side influent
SInPerm = DCMD_PackStream(SInPerm);
% 设定各级膜组件冷侧进料
SInPerms(1:NumStage) = SInPerm;
% 设定各级膜材料特性
Membranes(1:NumStage) = MembrProps;
% set properties for all TECs
load('TEC_Params.mat') % 载入已有的TEC计算参数
TECs(1:(NumStage+1)) = TEC_Params.TEC(3,1); % 注意按opt=0计算TEC的吸放热量
I1 = TECs(1).Current;
I2 = TECs(2).Current;

%% 计算SFMD性能
% 定义优化目标函数为单位能耗
opfun = @(x)interface_SFMD(x, NumStage, SInFeeds, SInPerms, Membranes, TECs);
% 初值
x0 = [W1 T1 I1 W2 T2 I2];
% 变量边界
lb = [3.5e-6, 40+273.15, 0, 3.5e-6, 5+273.15, 0];
ub = [2.1e-3, 80+273.15, 8, 2.1e-3, 30+273.15, 8];
% 求优
opts = optimoptions('fmincon','Display','iter');
x = fmincon(opfun,x0,[],[],[],[],lb,ub,[],opts);

%% 输出结果
[SEC,out1,out2] = opfun(x);
format short g
fprintf('Specific energy consumption of %d-stage DCMD is %.4e J/kg (%.4g kWh/kg).\n', NumStage, SEC, SEC*2.778e-7);
disp(out1)