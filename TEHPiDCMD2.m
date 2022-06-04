function [profile,SOUT] = TEHPiDCMD2(SINs,TECs,TEXs,membrane,flowPattern,opts)
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
nMesh = 21; % 包括两端点
% opt = optimset('Display','iter');
s2x0 = [mean([T1,T2]),W2];
s2lb = [273.2,W2];
s2ub = [T1,W2*10];
solOpts1 = optimoptions(@lsqnonlin,'Display','none');
f = @(x)CalcS2in(x)-[SINs(2).Temp,SINs(2).MassFlow];
s2x = lsqnonlin(f,s2x0,s2lb,s2ub,solOpts1);

% 输出结果
[~,profile,SOUT] = CalcS2in(s2x);

    % 求解膜组件逆流操作渗透侧出口温度和流率
    function [fval,profile,SOUT] = CalcS2in(x)
        s1in = SINs(1);
        s2out = SINs(2);
        s2out.MassFlow = x(2); 
        s2out.Temp = x(1);
        s2out.Enthalpy = s2out.MassFlow*s2out.SpecHeat*s2out.Temp;
        % 初值
        TSH = T1+1;
        TMF = mean([T1,mean([T1,T2])]);
        TMP = mean([T2,mean([T1,T2])]);
        TSC = TEXs(2);
        x0 = reshape([TSH;TMF;TMP;TSC]*ones(1,nMesh),[],1);
        % 求解温度侧形
        lb = ones(size(x0))*273.16;
        ub = ones(size(x0))*(273.15+98);
        fun = @(x)Eqns(x,TEXs,TECs,[s1in,s2out],membrane);
        % solOpts = optimoptions(@lsqnonlin,'Display','iter','Algorithm','trust-region-reflective');
        % [Ts,~,residual,exitflag] = lsqnonlin(fun,x0,lb,ub,solOpts);
        solOpts = optimoptions(@lsqnonlin,'Display','none','Algorithm','levenberg-marquardt');
        [Ts,~,residual,exitflag] = lsqnonlin(fun,x0,[],[],solOpts);
        [~,QTEC,QM,SM,S1,S2] = Eqns(Ts,TEXs,TECs,[s1in,s2out],membrane);
        fval(1) = S2(end).Temp;
        fval(2) = S2(end).MassFlow;
        if nargout > 1
            profile.T = reshape(Ts,[4,nMesh]);
            profile.T(:,end) = []; % 删除壁温最后一个无效节点，见Eqns()中iMesh循环
            profile.T = [profile.T(:,1),profile.T(:,1:end-1)+0.5*diff(profile.T,1,2)];
            S1T = [S1.Temp]; S1T(end) = [];
            S2T = [S2.Temp]; S2T(end) = [];
%             S1T = interp1(linspace(0,1,nMesh),[S1.Temp],linspace(0,1,nMesh-1));
%             S2T = interp1(linspace(0,1,nMesh),[S2.Temp],linspace(0,1,nMesh-1));
            profile.T = [profile.T;S1T;S2T];
            profile.T = profile.T([1,5,2,3,6,4],:);
            profile.QTEC = QTEC;
            profile.QM = QM;
            profile.SM = SM;
            profile.S1 = S1;
            profile.S2 = S2;
%             profile.Remarks = sprintf('DCMD膜组件流动模式：%s；exitflag（%d）',flowPattern,exitflag);
            profile.Remarks = sprintf('DCMD膜组件流动模式：%s',flowPattern);
            SOUT(1) = S1(end);
            SOUT(2) = S2(1);
        end
    end

    % 能量平衡方程组
    function [dQ,QTEC,QM,SM,S1,S2] = Eqns(Ts,TEXs,TECs,SINs,membrane)
        % 参数值
        AT(1:2) = 0.006*0.040; % 膜组件料液侧(1)和渗透侧(2)切向面积
        AN(1:2) = 0.040*0.040/(nMesh-1); % 膜组件料液侧(1)和渗透侧(2)法向面积
        x = reshape(Ts,[4,nMesh]);
        y = zeros(size(x));
        QTEC = cell(1,nMesh-1);
        QM = zeros(1,nMesh-1);
        SM(1:nMesh-1) = SINs(2);
        S1(1:nMesh) = SINs(1);
        S2(1:nMesh) = SINs(2);
        JM = zeros(1,nMesh);
        JH = zeros(1,nMesh);
        [~,KF] = DCMD_TM(SINs(1),0);
        [~,KP] = DCMD_TM(SINs(2),0);
        for iMesh = 1:nMesh-1
            QTEC1 = TE_Heat(x(1,iMesh),TEXs(1),TECs(1),opts(1),opts(2));
            QTEC2 = TE_Heat(TEXs(2),x(4,iMesh),TECs(2),opts(1),opts(2));
            QNF = QTEC1(1)/TECs(1).HTArea*AN(1);
            QSP = QTEC2(2)/TECs(2).HTArea*AN(2);
            QTEC{iMesh} = [QTEC1/TECs(1).HTArea*AN(1);QTEC2/TECs(2).HTArea*AN(2)];
            [JH(iMesh),JM(iMesh)] = JHM(SINs,x(2,iMesh),x(3,iMesh),membrane);
            QSF = JH(iMesh)*AN(1);
            QWF = S1(iMesh).Enthalpy;
            QEF = QWF+QNF-QSF;
            S1(iMesh+1).Enthalpy = QEF;
            S1(iMesh+1).MassFlow = S1(iMesh).MassFlow+JM(iMesh)*AN(1);
            S1(iMesh+1).Temp = S1(iMesh+1).Enthalpy/S1(iMesh+1).MassFlow/S1(iMesh+1).SpecHeat;
            TF = mean([S1(iMesh).Temp,S1(iMesh+1).Temp]);
            y(1,iMesh) = KF*(x(1,iMesh)-TF)-QNF;
            y(2,iMesh) = KF*(TF-x(2,iMesh))-QSF;
            QNP = JH(iMesh)*AN(2);
            QWP = S2(iMesh).Enthalpy;
            QEP = QWP+QSP-QNP;
            S2(iMesh+1).Enthalpy = QEP;
            S2(iMesh+1).MassFlow = S2(iMesh).MassFlow-JM(iMesh)*AN(2);
            S2(iMesh+1).Temp = S2(iMesh+1).Enthalpy/S2(iMesh+1).MassFlow/S2(iMesh+1).SpecHeat;
            TP = mean([S2(iMesh).Temp,S2(iMesh+1).Temp]);
            y(3,iMesh) = KP*(x(3,iMesh)-TP)-QNP;
            y(4,iMesh) = KP*(TP-x(4,iMesh))-QSP;            
            QM(iMesh) = mean([JH(iMesh)*AN(1),JH(iMesh)*AN(2)]);
            SM(iMesh).Temp = mean([x(2,iMesh),x(3,iMesh)]);
            SM(iMesh).MassFlow = JM(iMesh)*mean([AN(1),AN(2)]);
            SM(iMesh) = DCMD_PackStream(SM(iMesh),DuctGeom);
        end
        dQ = reshape(y,size(Ts));
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
%         % 计算壁面到主流的边界层内传热系数k及膜面边界层温差
%         switch flag
%             case 'feed-side'
%                 [~,k] = DCMD_TM(sIn(1),JH_BL);
%                 dT = sIn(1).Temp-TMF;
%             case 'perm-side'
%                 [~,k] = DCMD_TM(sIn(2),JH_BL);
%                 dT = TMP-sIn(2).Temp;
%         end
%         % 计算壁面到主流的边界层内传热通量
%         JH = k*dT;
    end
end