
function [outTab,profile] = SimDCMD(W1,T1,W2,T2)
outTab = table;
% 调用公用变量定义，其中包括DuctGeom（流道几何尺寸）、Stream（物料定义）、MembrProps（膜材料性质）
[DuctGeom,Stream,MembrProps] = InitStruct();
% 设定膜组件的热、冷侧进口温度和流率
if nargin == 0
    T1 = 323.15; T2 = 288.15; % [K]
    W1 = 1.217e-5; W2 = 1.217e-3; % [kg/s]
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
opts = [0,0]; TECs(1:(NumStage+1)) = TEC_Params.TEC(1,1); % 相当于未集成半导体热泵的DCMD膜组件


%% 计算集成热泵DCMD膜组件中的温度分布
[profile,~] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"countercurrent",opts);
opStr = 'cooling';

%% DCMD系统单位能耗
% 计算稳态操作时料液放热量Q(1)和渗透液吸热量Q(2)
[Q,WP,QM,TP1,TP2] = CalcHeat(profile);
outTab.WP = WP;
outTab.QM = QM;
% 加热器功耗
E(1) = Q(1);
outTab.Q1 = Q(1);
outTab.E1 = E(1);
outTab.Q2 = Q(2);
% 计算膜组件外置半导体制冷功耗
% opts = [0,1]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(18,1);
opts = [1,0]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(14,1);
% opts = [0,0]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(4,1);
[E(2),QTEC,~] = CalcTECPower(opStr,Q(2),TEXs(1),mean([TP1,TP2]),exTECs(1),opts);
outTab.QTEC = QTEC;
outTab.E2 = E(2);
% 计算系统总能耗
SEC = sum(E)/WP/3600/1000; % [kWh/kg]
outTab.SEC = SEC;

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





end