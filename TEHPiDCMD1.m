function [profile,SOUTs] = TEHPiDCMD1(SINs,TECs,TEXs,membrane,flowPattern,opts)
% 程序功能：
% 给定集成半导体热泵的DCMD膜组件料液和渗透侧进料流体状态（即SINs(1)和SINs(2)）；
% 膜组件两侧边界的半导体热泵操作条件TECs、设定参数opts及其冷热侧温度TEXs；
% 膜特性membrane；并或逆流操作模式flowPattern
% 计算膜组件内侧型profile；料液和渗透侧出料流体状态（即SOUTs(1)和SOUTs(2)）
%
% 计算说明：
% 膜组件中膜为水平放置，膜上侧为料液侧（代号为1）；下侧为渗透侧（代号为2）
% 料液侧流体固定从西向东流动；渗透侧流体从东向西流动，即逆流操作
% 沿流向将膜组件空间分为nMesh个网格，由此构成nMesh-1个控制体（CV）
% 依次从西到东对每个CV进行物热衡算，通过控制体进料焓计算出料焓，同时确定温度分布
%
% 调试用输入参数
if nargin == 0
    [DuctGeom,s1,membrane] = InitStruct();
    flowPattern = 'countercurrent';
    T1 = 273.15+50; T2 = 273.15+45; % [K]
    W1 = 1.217e-2; W2 = 6.389e-3; % [kg/s]
    s1.Temp = T1;
    s1.MassFlow = W1;
    SINs(1) = DCMD_PackStream(s1); % 膜组件料液侧进料
    s2 = s1;
    s2.Temp = T2;
    s2.MassFlow = W2;
    SINs(2) = DCMD_PackStream(s2); % 膜组件渗透侧进料
    iTECs = [23,21];
    T0 = 298.15;
    TEXs = [T0,T0];
    % set properties for all TECs
    load('TEC_Params.mat','TEC_Params') % 载入已有的TEC计算参数
    TECs(1) = TEC_Params.TEC(iTECs(1));
    % TEC计算参数
    opts = [TEC_Params.opt1(iTECs(1)),TEC_Params.opt2(iTECs(1))];
    TECs(2) = TEC_Params.TEC(iTECs(2));
else
    [DuctGeom,~,~] = InitStruct();
    T1 = SINs(1).Temp; T2 = SINs(2).Temp;
    W1 = SINs(1).MassFlow; W2 = SINs(2).MassFlow;
end
% 参数
sol1Opt = optimoptions(@lsqnonlin,'Display','none',...
                                 'Algorithm','trust-region-reflective');
sol2Opt = optimoptions(@lsqnonlin,'Display','none',...
                                 'Algorithm','levenberg-marquardt');
% 计算初值
nMesh = 2; % 包括两端点
x0 = [SINs(2).MassFlow,mean([SINs.Temp])];
% 注意用levenberg-marquardt算法会造成x(1)负值而出现复数解
lb = [0,273.15+5];
ub = [W2*10,273.15+95];
[x0,~,~,exitflag] = lsqnonlin(@CalcS2out,x0,lb,ub,sol1Opt);
if exitflag == 0
    prompt1 = sprintf('【注意】exitflag = 0');
else
    prompt1 = sprintf('exitflag = %d',exitflag);
    % 计算集成半导体热泵膜组件中的温度侧形
    nMesh = 11; % 似乎增大网格数会增加能量平衡偏差，例如设nMesh=21时的计算结果用ChkProfile()检查发现项目2、3偏差大于设定值
    QTEC = cell(1,nMesh-1);
    Ts = zeros(6,nMesh-1);
    S1(1:nMesh) = SINs(1);
    S2(1:nMesh) = SINs(2);
    SM(1:nMesh-1) = SINs(1);
    QM = zeros(1,nMesh-1);
%     x = lsqnonlin(@CalcS2out,x0,[],[],sol2Opt); % 逆流操作计算满足渗透侧进口条件的出口流率及温度; 
    lb = [0,273.15+5];
    ub = [W2*10,273.15+95];
    [x,~,res] = lsqnonlin(@CalcS2out,x0,lb,ub,sol1Opt);    
end
prompt = sprintf('%s；DCMD膜组件流动模式：%s',prompt1,flowPattern);
profile.T = Ts;
profile.QTEC = QTEC;
profile.QM = QM;
profile.SM = SM;
profile.S1 = S1;
profile.S2 = S2;
profile.Remarks = prompt;
SOUTs(1) = S1(end);
SOUTs(2) = S2(1);

    function y = CalcS2out(x) % 给定集成半导体热泵的DCMD膜组件渗透侧出口流率和温度
        nCV = nMesh-1;
        SWs(:,1) = SINs;
        SWs(2,1).MassFlow = x(1);
        SWs(2,1).Temp = x(2);
        SWs(2,1) = DCMD_PackStream(SWs(2,1),DuctGeom);
        SEs(:,nCV) = SINs;
        Ts = zeros(6,(nMesh-1));
        for iCV = 1:nCV
            [Ts(:,iCV),SEs(:,iCV),SM(iCV),QM(iCV),QTEC{iCV}] = HBCV(SWs(:,iCV));
            SWs(:,iCV+1) = SEs(:,iCV);
        end
        y(1) = SEs(2,end).MassFlow-SINs(2).MassFlow;
        y(2) = SEs(2,end).Temp-SINs(2).Temp;
        S1 = SWs(1,:);
        S2 = SWs(2,:);
    end

    function [Ts,SOUTs,SM,QM,QTEC] = HBCV(SINs)
        % 求解CV中满足进料条件的料液侧壁温、料液侧膜面温度、渗透侧膜面温度和渗透侧壁温Ts
        Ts0 = [T1 mean([T1,mean([T1,T2])]) mean([T2,mean([T1,T2])]) T2];
        f1 = @(x)CalcHB(x,SINs);
        [Ts,~,residual] = lsqnonlin(f1,Ts0,[],[],sol2Opt);
%         lb = 273.2*ones(size(Ts0));
%         ub = 372.2*ones(size(Ts0));
%         Ts = lsqnonlin(f1,Ts0,lb,ub,sol1Opt);
        [~,SOUTs,Ts,SM,QM,QTEC] = CalcHB(Ts,SINs);
    end

    function [y,SEs,Ts,SM,QM,QTEC] = CalcHB(x,SWs)
        % 给定CV中温度x及CV西侧流体状态SWs计算能量偏差y
        TSH = x(1); TMF = x(2); TMP = x(3); TSC = x(4);
        QFW = SWs(1).Enthalpy; QPW = SWs(2).Enthalpy;
        SM = SWs(2);
        AN(1:2) = 0.040*0.040/(nMesh-1); % 膜组件料液侧(1)和渗透侧(2)法向面积
        [~,KF] = DCMD_TM(SWs(1),0);
        [~,KP] = DCMD_TM(SWs(2),0);
        QTEC1 = TE_Heat(TSH,TEXs(1),TECs(1),opts(1),opts(2))/TECs(1).HTArea*AN(1);
        QTEC2 = TE_Heat(TEXs(2),TSC,TECs(2),opts(1),opts(2))/TECs(1).HTArea*AN(1);
        QFN = QTEC1(1);
        QPS = QTEC2(2);
        TF = TSH-QFN/AN(1)/KF;
        TP = TSC+QPS/AN(2)/KP;
        QFS = KF*(TF-TMF)*AN(1);
        QPN = KP*(TMP-TP)*AN(2);
        [JH,JM] = JHM(SWs,TMF,TMP,membrane);
        SEs = SWs;
        SEs(1).MassFlow = SWs(1).MassFlow-JM*AN(1);
        SEs(2).MassFlow = SWs(2).MassFlow-JM*AN(2);
        SEs(1).Temp = TF;
        SEs(2).Temp = TP;
        SEs(1) = DCMD_PackStream(SEs(1),DuctGeom);
        SEs(2) = DCMD_PackStream(SEs(2),DuctGeom);
        QFE = SEs(1).Enthalpy;
        QPE = SEs(2).Enthalpy;
        y(1) = QFW+QFN-QFE-QFS; % 料液侧CV进出能量偏差
        y(2) = QFS-JH*AN(1); % 料液侧膜面传热量与从CV中流向膜面热量的偏差
        y(3) = QPN-JH*AN(2); % 渗透侧膜面传热量与从膜面向CV中流出热量的偏差
        y(4) = QPE+QPN-QPW-QPS; % 渗透侧CV进出能量偏差
        Ts = [TSH,TF,TMF,TMP,TP,TSC]; % CV中的温度分布
        SM.MassFlow = JM*AN(2); % CV中跨膜流量
        SM.Temp = TMP; 
        SM = DCMD_PackStream(SM,DuctGeom);
        QM = QPN; % CV中跨膜传热量
        QTEC = [QTEC1;QTEC2];
    end
    
    % DCMD膜组件料液侧主流向膜面边界层的传热通量
    function [JH,JM] = JHM(sIn,TMF,TMP,membrane)
        MF = sIn(1).MassFraction;
        C  = membrane.MDCoefficient;
        K  = membrane.ThermConductivity;
        d  = membrane.Thickness;
        % 计算膜面蒸汽压
        PSH = arrayfun(@(x)DCMD_SatVapPressure(x,MF),TMF);
        PSC = arrayfun(@(x)DCMD_SatVapPressure(x),TMP);
        JM  = C*(PSH-PSC);
        % get temperature-dependent latent heat
        dHv = arrayfun(@(x)DCMD_LatentHeat(x),(TMF+TMP)/2);
        % 计算膜面传热通量
        JH = JM.*dHv+K/d*(TMF-TMP);
    end
end