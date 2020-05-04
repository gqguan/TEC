%% 回归DCMD过程计算参数（膜蒸馏系数）
%
% 功能说明
% （1）载入实验数据
% （2）回归膜蒸馏系数
%  调用DCMD求解器（将TEC设为近似的绝热壁面）计算相应实验设定条件下的膜组件出口温度和产水率
%
% by Dr. Guan Guoqiang @ SCUT on 2020-05-02
%

%% 初始化
% 载入数据结构
CommonDef
% 载入DCMD实验数据
import_DCMD_ExpData
ExpData = ExpData_raw(1:11,[9 10 11 12 13 15 17 19 21 23]);
% 将实验数据中的温度单位转换为K
ExpData.TF_IN = ExpData.TF_IN+273.15;
ExpData.TP_IN = ExpData.TP_IN+273.15;
ExpData.TF_OUT = ExpData.TF_OUT+273.15;
ExpData.TP_OUT = ExpData.TP_OUT+273.15;
% 将实验数据中的产水率单位转换为kg/s
ExpData.WP = ExpData.WP/1000;
% 设定集成TEC多级DCMD系统的级数
NumStage = 1;
% set properties for all TECs
load('TEC_Params.mat') % 载入已有的TEC计算参数
% 设定TEC输入电能为0、TEC导热系数为很小的值以此将TEC近似为绝热壁面
TEC_Params.TEC(3,1).Voltage = 0;
TEC_Params.TEC(3,1).Current = 0;
TEC_Params.TEC(3,1).ThermConductance = 1e-8;
TECs(1:(NumStage+1)) = TEC_Params.TEC(3,1); % 注意按opt=0计算TEC的吸放热量
% 设定冷热侧TEC的环境温度分别均为298.15K
TEXs = [298.15; 298.15];
% 求解结果变量
TSH = zeros(height(ExpData),1);
TH = zeros(size(TSH));
TMH = zeros(size(TSH));
TMC = zeros(size(TSH));
TC = zeros(size(TSH));
TSC = zeros(size(TSH));
WP_Sim = zeros(size(TSH));

%% 顺次求解各实验输入条件下的模组件出口温度和产水量
for i = 1:height(ExpData)
    % 膜组件热侧进料初始化
    T1 = ExpData.TF_IN(i); T2 = ExpData.TP_IN(i);
    SInFeed = Stream;
    SInFeed.Temp = T1;
    SInFeed.MassFlow = ExpData.QF_IN(i)*1e-3; % 注意这里将料液密度近似为1e-3 kg/mL
    % calculate the rest properties of feed-side influent
    SInFeed = DCMD_PackStream(SInFeed);
    % 设定各级膜组件热侧进料
    SInFeeds(1:NumStage) = SInFeed;
    % 膜组件冷侧进料初始化
    SInPerm = Stream;
    SInPerm.Temp = T2;
    SInPerm.MassFlow = ExpData.QP_IN(i)*1e-3; % 注意这里将渗透液密度近似为1e-3 kg/mL
    % calculate the rest properties of permeate-side influent
    SInPerm = DCMD_PackStream(SInPerm);
    % 设定各级膜组件冷侧进料
    SInPerms(1:NumStage) = SInPerm;
    % 设定各级膜材料特性
    Membranes(1:NumStage) = MembrProps;
    % 设定各级膜组件中的内部温度分布（热侧TEC壁面温度、热侧主体温度、热侧膜面温度、冷侧膜面温度、冷侧主体温度、冷侧TEC壁面温度）
    for j=1:NumStage
        T0((1+(j-1)*6):6*j) = [T1+1; T1; T1-1; T2+1; T2; T2-1];
    end
    % 求解满足能量平衡的模组件温度（热、冷侧的壁面温度、主体温度和膜面温度）
    opts = optimoptions('fsolve', 'Display', 'None', 'MaxFunEvals', 15000, 'MaxIter', 1000);
    fun = @(T)DCMD_EqSys(T, TEXs, TECs, SInFeeds, SInPerms, Membranes);
    [T, fvals, exitflag] = fsolve(fun, T0, opts);
    % 计算DCMD膜组件中的传热、传质量
    [~, Q, QM, SM, SOutFeeds, SOutPerms] = fun(T);
    % 整理温度向量
    TOut = reshape(T, [6,length(SM)]); % Temperature profiles of each stage
    TSH(i) = TOut(1,:)';
    TH(i)  = TOut(2,:)';
    TMH(i) = TOut(3,:)';
    TMC(i) = TOut(4,:)';
    TC(i)  = TOut(5,:)';
    TSC(i) = TOut(6,:)';
    % 计算产水量
    WP_Sim(i) = sum([SM.MassFlow]);
end

%% Output results
format short g
tabout = [ExpData(:,{'TF_OUT', 'TP_OUT', 'WP'}),table(TH, TC, WP_Sim)];
disp(tabout)
fprintf('RMS of relative errors in TF, TP and WP are %f, %f and %f, respectively\n', ...
rms((tabout.TF_OUT-tabout.TH)./tabout.TF_OUT), ...
rms((tabout.TP_OUT-tabout.TC)./tabout.TP_OUT), ...
rms((tabout.WP-tabout.WP_Sim)./tabout.WP))
% 画出实验与模拟的偏差图（散点横坐标为实验结果而纵坐标为模拟值，散点与对角线的距离表征实验与模拟的偏差，偏差越小距离越近）
figure
WP_max = max([tabout.WP;tabout.WP_Sim]);
WP_min = min([tabout.WP;tabout.WP_Sim]);
plot(tabout.WP, tabout.WP_Sim, 'o', [WP_min,WP_max], [WP_min,WP_max], '--r')
ax1 = gca;
ax1.XLim = [WP_min,WP_max];
ax1.YLim = [WP_min,WP_max];