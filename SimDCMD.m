
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
[~,Stream,MembrProps] = InitStruct();
% 设定膜组件的热、冷侧进口温度和流率
if ~exist('refluxRatio','var')
    refluxRatio = inf;
end
if ~exist('config','var')
    config = 'classical';
end
if nargin == 0
    T1 = 273.15+50; T2 = 273.15+45; % [K]
    W1 = 1.217e-2; W2 = 6.389e-3; % [kg/s]
    refluxRatio = inf; % 全回流
    config = 'feedTEHP'; 
end
% 设定集成TEC多级SFMD系统的级数
NumStage = 1;
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
        opts = [TEC_Params.opt1(iTEC1),TEC_Params.opt2(iTEC1)];
    case {'feedTEHP'} % DCMD膜组件料液侧集成半导体热泵
        iTEC1 = 9;
        TECs(1) = TEC_Params.TEC(iTEC1);
        opts = [TEC_Params.opt1(iTEC1),TEC_Params.opt2(iTEC1)];
        iTEC2 = 1; % 近似绝热
        TECs(2) = TEC_Params.TEC(iTEC2); 
        TEXs(1) = T2; % 即TECs(1)的冷侧温度初设为渗透液进膜组件温度
    otherwise
        error('SimDCMD()中输入参数config无法识别！')
end


%% 计算集成热泵DCMD膜组件中的温度分布
dTF2 = 1;
while abs(dTF2)>1e-8
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
            iTEC1 = 9;
            exTECs(1) = TEC_Params.TEC(iTEC1);
            opts = [TEC_Params.opt1(iTEC1),TEC_Params.opt2(iTEC1)];
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
            dTF2 = 0;
        case('extTEHP') % 计算稳态操作全回流时料液放热量Q(1)和渗透液吸热量Q(2)
            % 计算零回流时的料液加热所需热量Q(1)及渗透液冷却所需冷量Q(2)
            refluxRatio = 0;
            [Q,QM,WF,WP,TP1,TP2,TF1,TF2] = CalcHeat(profile,refluxRatio);
            % 计算膜组件外置半导体制冷功耗
            iTEC1 = 9;
            exTECs(1) = TEC_Params.TEC(iTEC1);
            opts = [TEC_Params.opt1(iTEC1),TEC_Params.opt2(iTEC1)];
            % 膜组件外集成TEC的冷热侧温度
            WR1 = refluxRatio*(WF-WP);
            TF0 = (WR1*TF2+WF*T0)/W1;
            TH = mean([TF0,TF1]);
            TC = mean([TP1,TP2]);
            % 计算满足渗透液冷却所需冷量的TEC输入电功
            [E(2),QTEC,~,nTEC] = CalcTECPower(opStr,Q(2),TH,TC,exTECs(1),opts);
            % TEC向料液侧放热量
            Q1 = E(2)+QTEC;
            if Q1 > Q(1) 
                warning('TEC放热量%.4g[W]大于料液加热所需最大热量%.4g[W]，建议提高W1！',Q1,Q(1))
            end
            % 计算料液侧加热量为Q1时的回流比
            [RR,~,WF,WP,~,~] = CalcReflux(profile,Q1);
            outTab.RR = RR;
            % 计算给定回流比时的TEC热侧温度
            WR1 = RR*(WF-WP);
            TF0 = (WR1*TF2+WF*T0)/W1;
            TH = mean([TF0,TF1]);
            % 计算满足渗透液冷却所需冷量的TEC输入电功
            [E(2),QTEC,~,nTEC] = CalcTECPower(opStr,Q(2),TH,TC,exTECs(1),opts);
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
            % 计算系统总能耗
            SEC = sum(E)/WP/3600/1000; % [kWh/kg]
            outTab.SEC = SEC;
            dTF2 = 0;
        case('feedTEHP')
            % 计算零回流时的料液加热所需热量
            refluxRatio = 0;
            [Q,QM,~,WP,TP1,TP2,TF1,TF2,dQ2] = CalcHeat(profile,refluxRatio,config);
            TH = mean([TF1,TF2]);
            TC = mean([TP1,TP2]);
            TEXs(1) = TC;
            % 按渗透液吸热量计算DCMD膜组件料液侧集成半导体热泵所需电功
            [TECs,profile1] = CalcTEHP(config,Q(2),sIn,TECs,TEXs,membrane,"countercurrent",opts);
            dTF2 = profile1.S1(end).Temp-TF2;
            iStart = strfind(profile.Remarks,'：');
            switch profile.Remarks(iStart+1:end)
                case('cocurrent')
                    dTP2 = profile1.S2(end).Temp-TP2;
                case('countercurrent')
                    dTP2 = profile1.S2(1).Temp-TP2;
            end
            if abs(dTF2) < 1e-8
                % 计算WF和R
                [RR,QM,WF,WP,~,~] = CalcReflux(profile1,Q(1));
                QTEC = sum(cellfun(@(x)x(1,2),profile1.QTEC));
                E(2) = sum(cellfun(@(x)x(1,1),profile1.QTEC))-QTEC;
                outTab.QTEC = QTEC;
                outTab.NTEC = 1;
                outTab.Q1 = Q(1);
                outTab.E1 = 0;
                outTab.Q2 = Q(2);
                outTab.E2 = E(2);
                outTab.RR = RR;
                % 计算系统总能耗
                SEC = sum(E)/WP/3600/1000; % [kWh/kg]
                outTab.SEC = SEC; 
                break
            else
                fprintf('dTF2 = %.4g[K]；dTP2 = %.4g[K]；dQ2 = %.4g[W]\n',dTF2,dTP2,dQ2);
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

end