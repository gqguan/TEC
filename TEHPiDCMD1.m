function [profile,SOUT] = TEHPiDCMD1(initTProfile,SIN,TEXs,membrane)
% 调试用输入参数
if nargin < 4
    [~,s1,membrane] = InitStruct();
    T1 = 273.15+50; T2 = 273.15+45; % [K]
    W1 = 1.217e-2; W2 = 6.389e-3; % [kg/s]
    s1.Temp = T1;
    s1.MassFlow = W1;
    SIN(1) = DCMD_PackStream(s1); % 膜组件料液侧进料
    s2 = s1;
    s2.Temp = T2;
    s2.MassFlow = W2;
    SIN(2) = DCMD_PackStream(s2); % 膜组件渗透侧进料
    iTECs = [20,1];
    T0 = 298.15;
    TEXs = [T0,T0];
    % set properties for all TECs
    load('TEC_Params.mat','TEC_Params') % 载入已有的TEC计算参数
    TECs(1) = TEC_Params.TEC(iTECs(1));
    % TEC计算参数
    opts = [TEC_Params.opt1(iTECs(1)),TEC_Params.opt2(iTECs(1))];
    TECs(2) = TEC_Params.TEC(iTECs(2)); 
end
% 参数值
nMesh = 21; % 包括两端点
AT(1:2) = 0.006*0.040; % 膜组件料液侧(1)和渗透侧(2)切向面积
AN(1:2) = 0.040*0.040/(nMesh-1); % 膜组件料液侧(1)和渗透侧(2)法向面积
% 物性和状态
% for i = 1:2
%     rho(i) = SIN(i).Density;
%     cp(i) = SIN(i).SpecHeat;
%     v(i) = SIN(i).Velocity;
% end
% 初值
if exist('initTProfile','var')
    x0 = reshape(initTProfile,[],1);
else
    TSH = T1+1;
    TF = T1;
    TP =T2;
    TMF = mean([T1,mean([T1,T2])]);
    TMP = mean([T2,mean([T1,T2])]);
    TSC = TEXs(2);
    x0 = reshape([TSH;TF;TMF;TMP;TP;TSC]*ones(1,nMesh),[],1);
end
% 求解温度侧形
lb = ones(size(x0))*273.16;
ub = ones(size(x0))*(273.15+98);
% solOpts = optimoptions(@lsqnonlin,'Display','iter','Algorithm','trust-region-reflective');
solOpts = optimoptions(@lsqnonlin,'Display','iter','Algorithm','levenberg-marquardt');
fun = @(x)Eqns(x,TEXs,TECs,SIN,membrane);
% [Ts,~,residual,exitflag] = lsqnonlin(fun,x0,lb,ub,solOpts);
[Ts,~,residual,exitflag] = lsqnonlin(fun,x0,[],[],solOpts);
% 输出结果
[~,QTEC,QM,SM,S1,S2] = Eqns(Ts,TEXs,TECs,SIN,membrane);
profile.T = reshape(Ts,[6,nMesh]);
profile.QTEC = QTEC;
profile.QM = QM;
profile.SM = SM;
profile.S1 = S1;
profile.S2 = S2;
profile.Remarks = sprintf('exitflag = %d',exitflag);
SOUT(1) = S1(end);
SOUT(2) = S2(1);

    % 能量平衡方程组
    function [fval,QTEC,QM,SM,S1,S2] = Eqns(xs,TEXs,TECs,SINs,membrane)
        x = reshape(xs,[6,nMesh]);
        y = zeros(size(x));
        QTEC = cell(1,nMesh);
        QM = zeros(1,nMesh);
        SM(1:nMesh) = SINs(2);
        S1(1:nMesh) = SINs(1);
        S2(1:nMesh) = SINs(2);
        JM = zeros(1,nMesh);
        a2 = zeros(1,nMesh);
        a5 = zeros(1,nMesh);
        [~,KF] = DCMD_TM(SINs(1),0);
        [~,KP] = DCMD_TM(SINs(2),0);
        x(2,1) = SINs(1).Temp;
        x(5,nMesh) = SINs(2).Temp;
        for iMesh = 1:nMesh
            QTEC1 = TE_Heat(x(1,iMesh),TEXs(1),TECs(1),opts(1),opts(2));
            QTEC2 = TE_Heat(TEXs(2),x(6,iMesh),TECs(2),opts(1),opts(2));
            JHWF = QTEC1(1)/TECs(1).HTArea;
            JHWP = QTEC2(2)/TECs(2).HTArea;
            QTEC{iMesh} = [QTEC1/TECs(1).HTArea*AN(1);QTEC2/TECs(2).HTArea*AN(2)];
            [JHMF,JM(iMesh)] = JHM(SINs,x(3,iMesh),x(4,iMesh),membrane);
            JHMP = JHM(SINs,x(3,iMesh),x(4,iMesh),membrane);
            QM(iMesh) = mean([JHMF*AN(1),JHMP*AN(2)]);
            y(3,iMesh) = KF*(x(2,iMesh)-x(3,iMesh))-JHMF;
            y(4,iMesh) = KP*(x(4,iMesh)-x(5,iMesh))-JHMP;
            a2(iMesh) = S1(iMesh).MassFlow*S1(iMesh).SpecHeat;
            a5(iMesh) = S2(iMesh).MassFlow*S2(iMesh).SpecHeat;
            switch iMesh
                case 1
                    SM(iMesh).MassFlow = 0.5*JM(iMesh)*AN(1);
                    S1(iMesh+1).MassFlow = S1(iMesh).MassFlow-SM(iMesh).MassFlow;
                    S2(iMesh+1).MassFlow = S2(iMesh).MassFlow-SM(iMesh).MassFlow;
                    a2(iMesh+1) = S1(iMesh+1).MassFlow*S1(iMesh+1).SpecHeat;
                    a5(iMesh+1) = S2(iMesh+1).MassFlow*S2(iMesh+1).SpecHeat;
                    y(1,iMesh) = KF*(x(1,iMesh)-x(2,iMesh))-JHWF;
                    y(2,iMesh) = KP*(x(5,iMesh)-x(6,iMesh))-JHWP;
                    y(5,iMesh) = (0.5*JHWF*AN(1)+a2(iMesh)*x(2,iMesh))-(a2(iMesh+1)*x(2,iMesh+1)+0.5*JHMF*AN(1));
                    y(6,iMesh) = (0.5*JHMP*AN(2)+a5(iMesh+1)*x(5,iMesh+1))-(a5(iMesh)*x(5,iMesh)+0.5*JHWP*AN(2));
                case nMesh
                    SM(iMesh).MassFlow = 0.5*JM(iMesh)*AN(1);
                    a2(iMesh-1) = S1(iMesh-1).MassFlow*S1(iMesh-1).SpecHeat;
                    a5(iMesh-1) = S2(iMesh-1).MassFlow*S2(iMesh-1).SpecHeat;
                    y(1,iMesh) = KF*(x(1,iMesh)-x(2,iMesh))-JHWF;
                    y(2,iMesh) = KP*(x(5,iMesh)-x(6,iMesh))-JHWP;
                    y(5,iMesh) = (0.5*JHWF*AN(1)+a2(iMesh)*x(2,iMesh))-(a2(iMesh-1)*x(2,iMesh-1)+0.5*JHMF*AN(1));
                    y(6,iMesh) = (0.5*JHMP*AN(2)+a5(iMesh)*x(5,iMesh))-(a5(iMesh-1)*x(5,iMesh-1)+0.5*JHWP*AN(2));
                otherwise
                    SM(iMesh).MassFlow = JM(iMesh)*AN(1);
                    S1(iMesh+1).MassFlow = S1(iMesh).MassFlow-SM(iMesh).MassFlow;
                    S2(iMesh+1).MassFlow = S2(iMesh).MassFlow-SM(iMesh).MassFlow;
                    a2(iMesh-1) = S1(iMesh-1).MassFlow*S1(iMesh-1).SpecHeat;
                    a2(iMesh+1) = S1(iMesh+1).MassFlow*S1(iMesh+1).SpecHeat;
                    a5(iMesh-1) = S2(iMesh-1).MassFlow*S2(iMesh-1).SpecHeat;
                    a5(iMesh+1) = S2(iMesh+1).MassFlow*S2(iMesh+1).SpecHeat;
                    y(1,iMesh) = KF*(x(1,iMesh)-x(2,iMesh))-JHWF;
                    y(2,iMesh) = KP*(x(5,iMesh)-x(6,iMesh))-JHWP;
                    y(5,iMesh) = (JHWF*AN(1)+0.5*a2(iMesh-1)*x(2,iMesh-1))-(0.5*a2(iMesh+1)*x(2,iMesh+1)+JHMF*AN(1));
                    y(6,iMesh) = (JHMP*AN(2)+0.5*a5(iMesh+1)*x(5,iMesh+1))-(0.5*a5(iMesh-1)*x(5,iMesh-1)+JHWP*AN(2));
            end
            SM(iMesh).Temp = mean([x(3,iMesh),x(4,iMesh)]);
        end
        fval = reshape(y,size(xs));
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