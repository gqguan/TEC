
function [outTab,profile] = SimDCMD(W1,T1,W2,T2,config,refluxRatio)
% 输入参数config说明：
% classical - 传统DCMD系统：外置料液加热和渗透液冷却单元，其中加热用电热，冷却用TEC
% extTEHP   - 在传统DCMD系统的基础上外置加热和冷却采用半导体热泵耦合，
%             采用料液侧部分回流解决TEC放热量大于吸热量的问题，故在该设定下无需输入参数refluxRation
% feedTEHP  - 在膜组件料液侧集成TEHP单元：TEHP热侧在膜组件料液流道中加热料液，而TEHP冷侧从渗透液吸热
outTab = table;
% WF = missing;
% WP = missing;
% QM = missing;
% 调用公用变量定义，其中包括DuctGeom（流道几何尺寸）、Stream（物料定义）、MembrProps（膜材料性质）
[DuctGeom,Stream,MembrProps] = InitStruct();
% 设定膜组件的热、冷侧进口温度和流率
if ~exist('refluxRatio','var')
    refluxRatio = inf;
end
if ~exist('config','var')
    config = 'classical';
end
if nargin == 0
    T1 = 318.15; T2 = 285.65; % [K]
    W1 = 1.217e-4*5; W2 = 6.146e-3; % [kg/s]
    refluxRatio = inf; % 全回流
    config = 'feedTEHP'; 
end
% 设定集成TEC多级SFMD系统的级数
NumStage = 1;
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
load('TEC_Params.mat','TEC_Params') % 载入已有的TEC计算参数
% 注意按opt1=0,opt2=1计算TEC的吸放热量
% opts = [0,1]; TECs(1:(NumStage+1)) = TEC_Params.TEC(18,1);
% opts = [1,0]; TECs(1:(NumStage+1)) = TEC_Params.TEC(14,1);
% opts = [0,0]; TECs(1:(NumStage+1)) = TEC_Params.TEC(4,1);
switch config
    case {'classical','extTEHP'} % 相当于未集成半导体热泵的DCMD膜组件
        opts = [0,0]; TECs(1:(NumStage+1)) = TEC_Params.TEC(1,1); % 近似绝热
    case {'feedTEHP'} % 集成半导体热泵的DCMD膜组件
        opts = [1,0]; TECs(1) = TEC_Params.TEC(15,1); % H27（RMSE=0.855）
        TECs(1).Current = 1.2;
        TECs(NumStage+1) = TEC_Params.TEC(2,1); % 近似绝热
    otherwise
        error('SimDCMD()中输入参数config无法识别！')
end


%% 计算集成热泵DCMD膜组件中的温度分布
relDiffQ = 1;
while abs(relDiffQ)>1e-8
    % 逆流操作（因为通常逆流操作单位能耗更低）
    [profile,~] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"countercurrent",opts);
    % % 并流操作
    % [profile,~] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"cocurrent",opts);
    opStr = 'cooling';

    %% DCMD系统单位能耗
    switch config 
        case('classical')
            % 计算稳态操作回流比为R时料液放热量Q(1)和渗透液吸热量Q(2)
            outTab.RR = refluxRatio;
            [Q,QM,WF,WP,TP1,TP2] = CalcHeat(profile,refluxRatio);
            % 计算膜组件外置半导体制冷功耗
            opts = [1,0]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(15,1);
            [E(2),QTEC,~,nTEC] = CalcTECPower(opStr,Q(2),TEXs(1),mean([TP1,TP2]),exTECs(1),opts);
            outTab.QTEC = QTEC;
            outTab.E2 = E(2);
            outTab.NTEC = nTEC;
            % 加热器功耗
            E(1) = Q(1);
            outTab.Q1 = Q(1);
            outTab.E1 = E(1);
            outTab.Q2 = Q(2);
            % 计算系统总能耗
            SEC = sum(E)/WP/3600/1000; % [kWh/kg]
            outTab.SEC = SEC;
            relDiffQ = 0;
        case('extTEHP') % 计算稳态操作全回流时料液放热量Q(1)和渗透液吸热量Q(2)
            % 计算零回流时的料液加热所需热量
            refluxRatio = 0;
            [Q,QM,~,WP,TP1,TP2] = CalcHeat(profile,refluxRatio);
            % 计算膜组件外置半导体制冷功耗
            % opts = [0,1]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(18,1);
            opts = [1,0]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(15,1);
            % opts = [0,0]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(4,1);
            [E(2),QTEC,~,nTEC] = CalcTECPower(opStr,Q(2),TEXs(1),mean([TP1,TP2]),exTECs(1),opts);
            outTab.QTEC = QTEC;
            outTab.E2 = E(2);
            outTab.NTEC = nTEC;
            % TEC向料液侧放热量
            Q1 = E(2)+QTEC;
            if Q1 > Q(1) 
                warning('TEC放热量%.4g[W]大于料液加热所需最大热量%.4g[W]，建议提高W1！',Q1,Q(1))
            end
            outTab.Q1 = Q1;
            outTab.E1 = 0; % 集成热泵时放热能耗为0
            outTab.Q2 = Q(2);
            % 计算料液侧加热量为Q1时的回流比
            [RR,~,WF,~,~,~] = CalcReflux(profile,Q1);
            outTab.RR = RR;
            % 计算系统总能耗
            SEC = sum(E)/WP/3600/1000; % [kWh/kg]
            outTab.SEC = SEC;
            relDiffQ = 0;
        case('feedTEHP')
            % 计算零回流时的料液加热所需热量
            refluxRatio = 0;
            [Q,QM,~,WP,TP1,TP2,TF1,TF2,relDiffQ] = CalcHeat(profile,refluxRatio);
            % 按渗透液吸热量计算TEC所需电功
            [E(2),QTEC,~,nTEC,TECs(1)] = CalcTECPower(opStr,Q(2),mean([TF1,TF2]),mean([TP1,TP2]),TECs(1),opts);
            if nTEC > 1
                warning('DCMD膜组件集成的TEC功率不满足当前指定的进料温度和流率条件')
            end
            if abs(relDiffQ)>1e-8 % 修正TEC输入电功
                msg = sprintf('TEC输入电流为%.4g[A]',TECs(1).Current);
                fprintf('%s：渗透液冷却所需冷量%.4g[W]与TEC吸热量相对偏差为%.4g%%\n',msg,Q(2),relDiffQ*100);
            else % 计算WF和R
                [RR,QM,WF,WP,~,~] = CalcReflux(profile,Q(1));
                outTab.QTEC = QTEC;
                outTab.nTEC = nTEC;
                outTab.Q1 = Q(1);
                outTab.E1 = 0;
                outTab.Q2 = Q(2);
                outTab.E2 = E(2);
                outTab.RR = RR;
                % 计算系统总能耗
                SEC = sum(E)/WP/3600/1000; % [kWh/kg]
                outTab.SEC = SEC;
            end
    end
end
outTab.WF = WF;
outTab.WP = WP;
outTab.QM = QM;
% 整理输出表格各列顺序
colNames = {'RR' 'WF' 'WP' 'QM' 'Q1' 'E1' 'Q2' 'QTEC' 'E2' 'NTEC' 'SEC'};
sortedColI = cellfun(@(x)find(strcmp(x,outTab.Properties.VariableNames)),colNames);
outTab = outTab(:,sortedColI);

%% 输出
% fprintf('DCMD系统在膜组件进口温度分别为%.2f[K]和%.2f[K]、流率分别为%.4g[kg/s]和%.4g[kg/s]\n',T1,T2,W1,W2)
% fprintf('DCMD系统稳定操作需要加热量%.4g[W]和制冷量%.4g[W]\n',Q(1),Q(2))
% fprintf('系统制冷采用TEC（编号%s）%s，计算热量为%.4g[W]，其偏差为%.4g[W]\n',TEC_Params.pid{14},opStr,QTEC,QDiff)
% fprintf('DCMD系统产水率为%.4g[kg/h]，加热消耗能量%.4g[W]，冷却消耗能量%.4g[W]，单位能耗为%.4g[kWh/kg]\n',WP*3600,E(1),E(2),SEC)
% DispResults(profile,DuctGeom,membrane);
% figObj2 = DispResults(profile2,DuctGeom,membrane);

function [DuctGeom,Stream,MembrProps] = InitStruct()
    % 定义模组件流道几何尺寸
    %  DuctGeom - geometric parameters of flowing channel in MD module
    %           .Length (real array(2)) length along the flowing direction in 
    %                                   both sides of MD module [m]
    %           .Height (real array(2)) height of rectanglarly wetted perimeter
    %                                   in both sides of MD module [m]
    %           .Width (real array(2))  width of rectanglarly wetted perimeter
    %                                   in both sides of MD module [m]
    DuctGeom = struct('Length', 0.04,  ...
                     'Height', 0.006, ...
                     'Width',  0.04 );
    % 定义物料数据结构
    Stream = struct('Temp', 295.15,  ...
                    'MassFlow', 0.015,  ...
                    'MassFraction', 0,  ...
                    'Density', 1e3,     ...
                    'Viscosity', 1e-3,  ...
                    'SpecHeat', 4.18e3, ...
                    'ThermCond', 0.6);
    % 定义膜材料性质
    %  MembrProps.TMH: hot-side temperature of membrane [K]
    %            .TMC: cold-side temperature of membrane [K]
    %            .Area: effective area of membrane [m2]
    %            .Thickness: thickness of membrane [m]
    %            .MDCoefficient: MD coefficient [kg/s-m2-Pa]
    %            .ThermConductivity: thermal conductivity of membrane
    MembrProps = struct('TMH', [], 'TMC', [], 'Area', 0.0016, ...
                        'Thickness', 1.5e-4, 'MDCoefficient', 3.2e-7, ...
                        'ThermConductivity', (0.18*0.3+0.025*0.7));
end







end