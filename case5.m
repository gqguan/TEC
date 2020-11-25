%% Calculate temperatures of stacked DCMD module with TEHPs
%  For the scenario of N-stage DCMD integrated with N+1 TEHPs
%
%  by Dr. GUAN, Guoqiang @ SCUT on 2019-08-30
%

%% Initialize
% 调用公用变量定义，其中包括DuctGeom（流道几何尺寸）、Stream（物料定义）、MembrProps（膜材料性质）
CommonDef
% 设定集成TEC多级DCMD系统的级数
NumStage = 1;
% 设定膜组件的热、冷侧进口温度
T1 = 323.15; T2 = 303.15;
% 膜组件热侧进料初始化
SInFeed = Stream;
SInFeed.Temp = T1;
SInFeed.MassFlow = 0.00001;
% calculate the rest properties of feed-side influent
SInFeed = DCMD_PackStream(SInFeed);
% 设定各级膜组件热侧进料
SInFeeds(1:NumStage) = SInFeed;
% 膜组件冷侧进料初始化
SInPerm = Stream;
SInPerm.Temp = T2;
SInPerm.MassFlow = 0.00001;
% calculate the rest properties of permeate-side influent
SInPerm = DCMD_PackStream(SInPerm);
% 设定各级膜组件冷侧进料
SInPerms(1:NumStage) = SInPerm;
% 设定各级膜材料特性
Membranes(1:NumStage) = MembrProps;
% set properties for all TECs
load('TEC_Params.mat') % 载入已有的TEC计算参数
TECs(1:(NumStage+1)) = TEC_Params.TEC(3,1); % 注意按opt=0计算TEC的吸放热量
% Debugging ======= 禁用TEC加热和冷却
for iTEC = 1:length(TECs)
    TECs(iTEC).Voltage = 0;
    TECs(iTEC).Current = 0;
end
% ============
% 设定各级膜组件中的内部温度分布（热侧TEC壁面温度、热侧主体温度、热侧膜面温度、冷侧膜面温度、冷侧主体温度、冷侧TEC壁面温度）
for i=1:NumStage
    T0((1+(i-1)*6):6*i) = [T1+1; T1; T1-1; T2+1; T2; T2-1];
end
% set room temperature as the environmental temperature of both heat source
% and sink
TEXs = [298.15; 298.15];

%% Solve temperatures
opts = optimoptions('fsolve', 'Display', 'iter', 'MaxFunEvals', 15000, 'MaxIter', 1000);
fun = @(T)DCMD_EqSys(T, TEXs, TECs, SInFeeds, SInPerms, Membranes);
[T, fvals, exitflag] = fsolve(fun, T0, opts);
[~, Q, QM, SM, SOutFeeds, SOutPerms] = fun(T);

%% Calculate the energy efficiency
% Energy consumption of TECs
EC = Q(:,1)-Q(:,2);
% Specific energy consumption
WP_Sum = sum([SM.MassFlow]);
SEC = sum(EC)/WP_Sum;
% Coefficient of performance
COP_H = Q(:,1)./EC; % heating COPs of TEHP
COP_C = Q(:,2)./EC; % refrigrating COPs of TEHP

%% Output results
format short g
fprintf('Specific energy consumption of %d-stage DCMD is %.4e W/kg.\n', NumStage, SEC);
% Performances of membrane separation in each stage
StageNames = cell(NumStage,1);
JM = zeros(NumStage,1);
JH = zeros(NumStage,1);
for i = 1:NumStage
    StageNames{i} = sprintf('Stage-%d', i);
    JM(i) = SM(i).MassFlow/Membranes(i).Area;
    JH(i) = QM(i)/Membranes(i).Area;
end
TOut = reshape(T, [6,length(SM)]); % Temperature profiles of each stage
TSH = TOut(1,:)';
TH  = TOut(2,:)';
TMH = TOut(3,:)';
TMC = TOut(4,:)';
TC  = TOut(5,:)';
TSC = TOut(6,:)';
Output_Stages = table(TSH, TH, TMH, TMC, TC, TSC, JH, JM, ...
                      'RowNames', StageNames);
disp(Output_Stages)
% Performances of TEHPs
TECNames = cell(NumStage+1,1);
for i = 1:(NumStage+1)
    TECNames{i} = sprintf('TEHP-%d', i);
end
TS_Hots(1:NumStage,1)  = TOut(1,:)';
TS_Hots(NumStage+1,1) = TEXs(2);
TS_Colds(1,1) = TEXs(1);
TS_Colds(2:NumStage+1,1) = TOut(6,:)';
Q1 = Q(:,1);
Q2 = Q(:,2);
Output_TEHP = table(COP_H, COP_C, TS_Hots, TS_Colds, Q1, Q2, EC, ...
                    'RowNames', TECNames);
disp(Output_TEHP)