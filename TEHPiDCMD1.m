function [T,residual,exitflag] = TEHPiDCMD1(SIN,iTECs,TEXs,membrane)
% 调试用输入参数
if nargin == 0
    [~,s1,membrane] = InitStruct();
    T1 = 273.15+50; T2 = 273.15+45; % [K]
    W1 = 1.217e-2; W2 = 6.389e-3; % [kg/s]
%     config = 'permTEHP'; 
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
for i = 1:2
    rho(i) = SIN(i).Density;
    cp(i) = SIN(i).SpecHeat;
    v(i) = SIN(i).Velocity;
end
% 初值
TSH = T1+1;
TF = T1;
TP =T2;
TMF = mean([T1,mean([T1,T2])]);
TMP = mean([T2,mean([T1,T2])]);
TSC = TEXs(2);
x0 = reshape([TSH;TF;TMF;TMP;TP;TSC]*ones(1,nMesh),[],1);
lb = ones(size(x0))*273.16;
ub = ones(size(x0))*(273.15+98);
solOpts = optimoptions(@lsqnonlin,'Display','iter');
fun = @(x)Eqns(x,TEXs,TECs,SIN,membrane);
[T,~,residual,exitflag] = lsqnonlin(fun,x0,lb,ub,solOpts);
T = reshape(T,[6,nMesh]);

    % 能量平衡方程组
    function fval = Eqns(xs,TEXs,TECs,SINs,membrane)
        x = reshape(xs,[6,nMesh]);
        y = zeros(size(x));
        JM = zeros(1,nMesh);
        a2 = zeros(1,nMesh);
        a5 = zeros(1,nMesh);
        [~,KF] = DCMD_TM(SINs(1),0);
        [~,KP] = DCMD_TM(SINs(2),0);
        for iMesh = 1:nMesh
            QTEC1 = TE_Heat(x(1,iMesh),TEXs(1),TECs(1),opts(1),opts(2));
            QTEC2 = TE_Heat(TEXs(2),x(6,iMesh),TECs(2),opts(1),opts(2));
            JHWF = QTEC1(1)/TECs(1).HTArea;
            JHWP = QTEC2(2)/TECs(2).HTArea;
            [JHMF,JM(iMesh)] = JHM(SINs,x(3,iMesh),x(4,iMesh),membrane);
            JHMP = JHM(SINs,x(3,iMesh),x(4,iMesh),membrane);
            y(3,iMesh) = KF*(x(2,iMesh)-x(3,iMesh))-JHMF;
            y(4,iMesh) = KP*(x(4,iMesh)-x(5,iMesh))-JHMP;
            massFlow = AT.*rho.*v;
            a2(iMesh) = massFlow(1)*cp(1);
            a5(iMesh) = massFlow(2)*cp(2);
            switch iMesh
                case 1
                    a2(iMesh+1) = (massFlow(1)-JM(iMesh)*AN(1))*cp(1);
                    a5(iMesh+1) = (massFlow(2)-JM(iMesh)*AN(2))*cp(2);
                    y(1,iMesh) = KF*(x(1,iMesh)-x(2,iMesh))-2*JHWF;
                    y(2,iMesh) = KP*(x(5,iMesh)-x(6,iMesh))-2*JHWP;
                    y(5,iMesh) = JHWF*AN(1)+a2(iMesh)*x(2,iMesh)-a2(iMesh+1)*x(2,iMesh+1)-JHMF*AN(1);
                    y(6,iMesh) = JHMP*AN(2)+a5(iMesh+1)*x(5,iMesh+1)-a5(iMesh)*x(5,iMesh)-JHWP*AN(2);
                case nMesh
                    a2(iMesh-1) = (massFlow(1)+JM(iMesh)*AN(1))*cp(1);
                    a5(iMesh-1) = (massFlow(2)+JM(iMesh)*AN(2))*cp(2);
                    y(1,iMesh) = KF*(x(1,iMesh)-x(2,iMesh))-2*JHWF;
                    y(2,iMesh) = KP*(x(5,iMesh)-x(6,iMesh))-2*JHWP;
                    y(5,iMesh) = 2*JHWF*AN(1)+a2(iMesh)*x(2,iMesh)-a2(iMesh-1)*x(2,iMesh-1)-JHMF*AN(1);
                    y(6,iMesh) = JHMP*AN(2)+a5(iMesh)*x(5,iMesh)-a5(iMesh-1)*x(5,iMesh-1)-2*JHWP*AN(2);
                otherwise
                    a2(iMesh-1) = (massFlow(1)+JM(iMesh)*AN(1))*cp(1);
                    a2(iMesh+1) = (massFlow(1)-JM(iMesh)*AN(1))*cp(1);
                    a5(iMesh-1) = (massFlow(2)+JM(iMesh)*AN(2))*cp(2);
                    a5(iMesh+1) = (massFlow(2)-JM(iMesh)*AN(2))*cp(2);
                    y(1,iMesh) = KF*(x(1,iMesh)-x(2,iMesh))-JHWF;
                    y(2,iMesh) = KP*(x(5,iMesh)-x(6,iMesh))-JHWP;
                    y(5,iMesh) = JHWF*AN(1)+0.5*a2(iMesh-1)*x(2,iMesh-1)-0.5*a2(iMesh+1)*x(2,iMesh+1)-JHMF*AN(1);
                    y(6,iMesh) = JHMP*AN(2)+0.5*a5(iMesh+1)*x(5,iMesh+1)-0.5*a5(iMesh-1)*x(5,iMesh-1)-JHWP*AN(2);
            end
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