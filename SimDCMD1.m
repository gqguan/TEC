
function [outTab,profile] = SimDCMD1(sn,W1,T1,W2,T2,config,refluxRatio)
% 输入参数config说明：
% classical - 传统DCMD系统：外置料液加热和渗透液冷却单元，其中加热用电热，冷却用TEC
% extTEHP   - 在传统DCMD系统的基础上外置加热和冷却采用半导体热泵耦合，
%             采用料液侧部分回流解决TEC放热量大于吸热量的问题，故在该设定下无需输入参数refluxRation
% feedTEHP  - 在膜组件料液侧集成TEHP单元：TEHP热侧在膜组件料液流道中加热料液，而TEHP冷侧从渗透液吸热
outTab = table;
WF = nan;
WP = nan;
QM = nan;
TF0 = nan;
strIU = {'Current','Voltage'};
strUnit = {'A','V'};
% 调用公用变量定义，其中包括DuctGeom（流道几何尺寸）、Stream（物料定义）、MembrProps（膜材料性质）
[~,Stream,MembrProps] = InitStruct();
% 设定膜组件的热、冷侧进口温度和流率
if ~exist('refluxRatio','var')
    refluxRatio = inf;
end
if ~exist('config','var')
    config = 'classical';
end
if nargin == 0
    sn = 'test';
    T1 = 330.65; T2 = 318.15; % [K]
    W1 = 1.217e-4*25; W2 = 1.217e-2; % [kg/s]
    refluxRatio = inf; % 全回流
    flowPattern = "countercurrent";
    config = 'extTEHP'; 
end

% 设定环境温度
TEXs = [298.15,298.15]; % [K]
T0 = TEXs(1);
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

switch config
    case {'classical','extTEHP'} % 相当于未集成半导体热泵的DCMD膜组件
        iTEC1 = 1;
        TECs(1:2) = TEC_Params.TEC(iTEC1);
        TECOpts = [TEC_Params.opt1(iTEC1),TEC_Params.opt2(iTEC1)];
        % 计算膜组件内温度分布
        profile = TEHPiDCMD(sIn,TECs,TEXs,membrane,flowPattern,TECOpts);
        iStart = strfind(profile.Remarks,'：');
        switch profile.Remarks(iStart+1:end)
            case('cocurrent')
                TP1 = profile.S2(1).Temp;
                TP2 = profile.S2(end).Temp;
                W2 = profile.S2(1).MassFlow;
            case('countercurrent')
                TP1 = profile.S2(end).Temp;
                TP2 = profile.S2(1).Temp;
                W2 = profile.S2(end).MassFlow;
            otherwise
                error('CalcHeat()输入参数profile字段Remarks中无有效的流型信息')
        end
        % DCMD系统渗透侧所需冷量（对给定的膜组件内温度分布，料液侧回流比不影响渗透侧冷量）
        [Q,QM,~,WP,TP1,TP2,~,~,~,~] = CalcHeat(profile,inf,config);
        % 计算膜组件外置半导体制冷功耗
        iTEC1 = 20;
        exTECs(1) = TEC_Params.TEC(iTEC1);
        opts = [TEC_Params.opt1(iTEC1),TEC_Params.opt2(iTEC1)];
        iTEC = 0;
    case {'feedTEHP'} % DCMD膜组件料液侧集成半导体热泵
        iTEC1 = 20;
        TECs(1) = TEC_Params.TEC(iTEC1);
        opts = [TEC_Params.opt1(iTEC1),TEC_Params.opt2(iTEC1)];
        iTEC2 = 1; % 近似绝热
        TECs(2) = TEC_Params.TEC(iTEC2); 
        TEXs(1) = T2; % 即TECs(1)的冷侧温度初设为渗透液进膜组件温度
        iTEC = 1; % 用于识别TEC(1)还是TEC(2)的传热量计算：这里为TEC(1)
    case {'permTEHP','permTEHP1','permTEHP2'} % DCMD膜组件渗透侧集成半导体热泵
        iTEC1 = 1;
        TECs(1) = TEC_Params.TEC(iTEC1);
        iTEC2 = 20;
        TECs(2) = TEC_Params.TEC(iTEC2);
        opts = [TEC_Params.opt1(iTEC2),TEC_Params.opt2(iTEC2)];
        TEXs(2) = T1; % TECs(2)的热侧温度初设为料液进膜组件温度
        iTEC = 2; % 用于识别TEC(1)还是TEC(2)的传热量计算：这里为TEC(2)
    otherwise
        error('SimDCMD()中输入参数config无法识别！')
end

% 边界条件
lb = [273.15+5,0.2];
ub = [273.15+90,15];
% 求解器参数
solOpt = optimoptions(@lsqnonlin,'Display','iter');
% 目标函数求解
switch config
    case 'extTEHP'
        x0 = [T1,1.5]; % 初值[TEC热侧平均温度，TEC操作电流或电压值]
        f = @(x)FunSys1(x,exTECs(1),opts,Q(2),mean([TP1,TP2]),profile,config);
    case {'feedTEHP','permTEHP','permTEHP'}
        x0 = [T2,1.5]; % 初值[TEC冷侧平均温度，TEC操作电流或电压值]
        f = @(x)FunSys(x,sIn,TECs,TEXs,membrane,flowPattern,opts,config);
end
[xsol,~,residual,exitflag] = lsqnonlin(f,x0,lb,ub,solOpt);

%% 输出

if exitflag == 1
    switch config
        case 'extTEHP'
            fprintf('TEC(%d)热侧平均温度为%.4g[K]、操作%s为%.4g[%s]\n',...
                iTEC,xsol(1),strIU{opts(2)+1},xsol(2),strUnit{opts(2)+1})
            fprintf('TEC(%d)热侧水箱热量衡算偏差和平均温度偏差分别为%.4g[W]和%.4g[K]\n',...
                iTEC,residual)
            Q2 = Q(2);
        case {'feedTEHP','permTEHP','permTEHP1'}
            fprintf('TEC(%d)冷侧平均温度为%.4g[K]、操作%s为%.4g[%s]\n',...
                iTEC,xsol(1),strIU{opts(2)+1},xsol(2),strUnit{opts(2)+1})
            fprintf('TEC(%d)冷侧水箱热量衡算偏差和平均温度偏差分别为%.4g[W]和%.4g[K]\n',...
                iTEC,residual)
            Q1 = sum(cellfun(@(x)x(iTEC,1),profile.QTEC));
            Q2 = sum(cellfun(@(x)x(iTEC,2),profile.QTEC));
            [RR,~,~,~,~,~] = CalcReflux(profile,Q1);
            [Q,QM,WF,WP,TP1,TP2,TF1,TF2,dQ2] = CalcHeat(profile,RR,config);
            % 计算系统总能耗
            E(1) = 0;
            E(2) = Q1-Q2;
    end
    SEC = sum(E)/WP/3600/1000; % [kWh/kg]
    % 输出变量
    outTab.RR = RR;
    outTab.WF = WF;
    outTab.WP = WP;
    outTab.QM = QM;
    outTab.Q1 = Q(1);
    outTab.E1 = E(1);
    outTab.Q2 = Q(2);
    outTab.E2 = E(2);
    outTab.QTEC = Q2;
    outTab.NTEC = 1;
    outTab.TF0 = nan;
    outTab.SEC = SEC;
    outTab.NOTE = {'ok'};
end
% 整理输出表格各列顺序
colNames = {'RR' 'WF' 'WP' 'QM' 'Q1' 'E1' 'Q2' 'E2' 'QTEC' 'NTEC' 'TF0' 'SEC' 'NOTE'};
sortedColI = cellfun(@(x)find(strcmp(x,outTab.Properties.VariableNames)),colNames);
outTab = outTab(:,sortedColI);

    function fval = FunSys1(x,TEC,opts,Q2,Tc,profile,cfg)
        Th = x(1);
        TEC.(strIU{opts(2)+1}) = x(2);
        Q = TE_Heat(Th,Tc,TEC,opts(1),opts(2));
        [RR,~,~,~,~,~] = CalcReflux(profile,Q(1));
        [Q,QM,WF,WP,TP1,TP2,TF1,TF2,~,TF0] = CalcHeat(profile,RR,cfg);
        [E(2),~,~,nTEC,newTEC] = CalcTECPower('cooling',Q(2),Th,Tc,TEC,opts);
        fval(1) = Q(2)-Q2;
        fval(2) = x(1)-mean([TF0,TF1]);
    end

function fval = FunSys(x,sIn,TECs,TEXs,membrane,flowPattern,TECOpts,cfg)
    fval = zeros(size(x));
    switch cfg
        case 'feedTEHP'
            TECs(1).(strIU{TECOpts(2)+1}) = x(2);
            TEXs(1) = x(1); 
            profile = TEHPiDCMD(sIn,TECs,TEXs,membrane,flowPattern,TECOpts);
            iStart = strfind(profile.Remarks,'：');
            switch profile.Remarks(iStart+1:end)
                case('cocurrent')
                    TP1 = profile.S2(1).Temp;
                    TP2 = profile.S2(end).Temp;
                    W2 = profile.S2(1).MassFlow;
                case('countercurrent')
                    TP1 = profile.S2(end).Temp;
                    TP2 = profile.S2(1).Temp;
                    W2 = profile.S2(end).MassFlow;
                otherwise
                    error('CalcHeat()输入参数profile字段Remarks中无有效的流型信息')
            end
            % 膜组件渗透侧进出口温度
            cp2 = mean([profile.S2.SpecHeat]);
            Q2TEC = sum(cellfun(@(x)x(1,2),profile.QTEC));
            fval(1) = Q2TEC-W2*cp2*(TP2-TP1);
            fval(2) = x(1)-mean([TP1,TP2]);
        case {'permTEHP','permTEHP1'}
            TECs(2).(strIU{TECOpts(2)+1}) = x(2);
            TEXs(2) = x(1);
            profile = TEHPiDCMD(sIn,TECs,TEXs,membrane,flowPattern,TECOpts);
            iStart = strfind(profile.Remarks,'：');
            switch profile.Remarks(iStart+1:end)
                case('cocurrent')
                    TP1 = profile.S2(1).Temp;
                    TP2 = profile.S2(end).Temp;
                    W2 = profile.S2(1).MassFlow;
                case('countercurrent')
                    TP1 = profile.S2(end).Temp;
                    TP2 = profile.S2(1).Temp;
                    W2 = profile.S2(end).MassFlow;
                otherwise
                    error('CalcHeat()输入参数profile字段Remarks中无有效的流型信息')
            end
            % 膜组件料液加热单元出口温度
            cp1 = mean([profile.S1.SpecHeat]);
            Q1TEC = sum(cellfun(@(x)x(2,1),profile.QTEC));
            [RR,~,~,~,~,~] = CalcReflux(profile,Q1TEC);
            [Q,QM,WF,WP,TP1,TP2,TF1,TF2,dQ1,TF0] = CalcHeat(profile,RR,config);
            fval(1) = Q1TEC-W1*cp1*(TF1-TF0);
            fval(2) = x(1)-mean([TF1,TF0]);
    end
end

end