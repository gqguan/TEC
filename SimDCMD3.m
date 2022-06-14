function [stream,QTEC,profile,exitflag] = SimDCMD3(Eset,W1,W2,T0,cfg)
% 模拟给定功耗和冷热侧循环流率的外置集成半导体热泵DCMD系统
    % 初始化
    strIU = {'Current','Voltage'};
    strUnit = {'A','V'};
    % 获取DuctGeom（流道几何尺寸）、Stream（物料定义）、MembrProps（膜材料性质）
    [~,Stream,MembrProps] = InitStruct();
    flowPattern = 'countercurrent';
    % 设定外置集成半导体热泵换热单元的进料[stmF0,stmP2]初值
    if nargin == 0
        % 设定系统操作参数
        Eset = 24; % TEC功耗
        W1 = 1.217e-2; W2 = 6.389e-3; % [kg/s]
        T0 = 298.15;
        cfg = 'classical';
    end
    % 设定半导体热泵参数
    iTEC1 = 22;
    load('TEC_Params.mat','TEC_Params') % 载入已有的TEC计算参数
    extTEC = TEC_Params.TEC(iTEC1);
    opts = [TEC_Params.opt1(iTEC1),TEC_Params.opt2(iTEC1)];
    % 膜组件设定
    TECs(1) = TEC_Params.TEC(21); % 绝热边界
    TECs(2) = TEC_Params.TEC(21); % 绝热边界
    TEXs(1) = T0; % 膜组件的环境温度
    TEXs(2) = T0;
    membrane = MembrProps; % 膜性质
    % 设初值
    T1 = 273.15+50; T2 = 273.15+45; % 冷热侧温度 [K]
    WF = 1.5e-5; % 进料量 [kg/s]
    % 系统进料（stream.F）初值，变量stream.F会在后续迭代中修正
    stream.F = Stream;
    stream.F.Temp = T0;
    stream.F.MassFlow = WF;
    stream.F = DCMD_PackStream(stream.F);

    opt = optimoptions(@lsqnonlin,'Display','iter',...
                                  'Algorithm','trust-region-reflective');
    switch cfg
        case('classical')
            % 初始化流股F1和R2
            stream.F1 = Stream; stream.R2 = Stream;
            stream.F1.MassFlow = W1; stream.R2.MassFlow = W2;
            switch strIU{opts(2)}
                case('Current')
                    initIorU = 3.5;
                    ubIorU = 15;
                case('Voltage')
                    initIorU = 6;
                    ubIorU = 15;
            end
            % 非线性最小二乘法求解[TF1,TR2,IorU]
            x0 = [T1,T2,initIorU];
            lb = [T0,273.15+5,0.1];
            ub = [273.15+95,T1,ubIorU];
            [x,~,residual,exitflag] = lsqnonlin(@classical,x0,lb,ub,opt);
            stmList = {"F","P","F0","F1","F2","P1","P2","R2"};
            fprintf('求解TF1、TP1和TEC操作%s分别为%.5gK、%.5gK和%.5g%s（exitflag=%d）；\n',strIU{opts(2)},x,strUnit{opts(2)},exitflag)
            fprintf('其中残差分别为%.5gK、%.5gK和%.5g%s\n',residual,strUnit{opts(2)})
        case('extTEHP')
            % 初始化流股F0和R2
            stream.F0 = Stream; stream.R2 = Stream;
            stream.F0.MassFlow = W1; stream.R2.MassFlow = W2;
            TF0 = mean([T1,T2]); TP2 = T2;
            % 非线性最小二乘法求解[TF0,TR2]
            x0 = [TF0,T2];
            lb = [T0,273.15+5];
            ub = [273.15+95,T1];
            [x,~,residual,exitflag] = lsqnonlin(@extTEHP1,x0,lb,ub,opt);
            stmList = {"F","FP","P","F0","F1","F2","R1","P1","P2","R2"};
            fprintf('求解TF0和T2分别为%.5gK和%.5gK（exitflag=%d）；其中残差分别为%.5gK和%.5gK\n',x,exitflag,residual)
    end
    
%     % 直接迭代法
%     for iter = 1:1000 % 最大迭代次数为1000
%         [TF0a,TP2a] = extTEHP(TF0,TP2); 
%         disp([TF0a,TP2a])
%         if IsConvergence([TF0a,TF0],[TP2a,TP2])
%             break
%         else
%             TF0 = TF0a;
%             TP2 = TP2a;
%         end
%     end
%     if iter == 1000
%         exitflag = 0;
%     else
%         exitflag = 1;
%     end
    % 输出
    outTab = table;
    for iStm = 1:length(stmList)
        Stream = stmList{iStm};
        tab1r = [table(Stream),struct2table(stream.(stmList{iStm}))];
        outTab = [outTab;tab1r];
    end
    disp(outTab)
    
    function y = classical(x)
        y = ones(size(x)); % 设初值
        stream.F1.Temp = x(1); stream.R2.Temp = x(2); extTEC.(strIU{opts(2)}) = x(3);
        stream.F1 = DCMD_PackStream(stream.F1);
        stream.R2 = DCMD_PackStream(stream.R2);
        [stream.P1,QTEC] = TEC(stream.R2,extTEC,opts);
        E2 = -diff(QTEC);
        [profile,SOUTs] = TEHPiDCMD1([stream.F1,stream.P1],TECs,TEXs,membrane,flowPattern,opts);
        stream.F2 = SOUTs(1);
        stream.P2 = SOUTs(2);
        y(2) = (stream.P2.Temp-x(2))/x(2);
        stream.P = SPLIT(stream.P2,stream.R2);
        stream.F = CalcF(stream.P);
        E1 = CalcE1(stream.F1,stream.F2,stream.F);
        stream.F0 = MIX(stream.F2,stream.F);
        stmF1 = Heater(stream.F0,E1);
        y(1) = (stmF1.Temp-x(1))/x(1);
        y(3) = (EDistributor(Eset,E1)-E2)/E2;
    end

    function E1 = CalcE1(stmF1,stmF2,stmF)
        E1 = stmF1.Enthalpy-stmF2.Enthalpy-stmF.Enthalpy;
    end

    function [sout,QTEC] = TEC(sin,extTEC,opts)
        sout = sin; % 设初值
        % 计算TEC冷却渗透液出料
        TC0 = sin.Temp;
        TC = fzero(@CalcTC,TC0);
        
        function dTC = CalcTC(TC)
            QTEC = TE_Heat(298.15,TC,extTEC,opts(1),opts(2));
            sout.Enthalpy = sin.Enthalpy-QTEC(2);
            sout.Temp = sout.Enthalpy/sout.MassFlow/sout.SpecHeat;
            dTC = sout.Temp-TC;
        end
    end

    function stmF = CalcF(stmP)
        stmF = stmP;
        stmF.Temp = T0;
        stmF = DCMD_PackStream(stmF);
    end

    function [TF0a,TP2a] = extTEHP(TF0,TP2)
        stream.F0.Temp = TF0; stream.R2.Temp = TP2;
        stream.F0 = DCMD_PackStream(stream.F0);
        stream.R2 = DCMD_PackStream(stream.R2);
        [SOUTs,QTEC] = HXiTEHP([stream.F0,stream.R2],extTEC,opts,Eset);
        stream.F1 = SOUTs(1);
        stream.P1 = SOUTs(2);
        SINs = SOUTs;
        [profile,SOUTs] = TEHPiDCMD1(SINs,TECs,TEXs,membrane,flowPattern,opts);
        stream.F2 = SOUTs(1);
        stream.P2 = SOUTs(2);
        TP2a = stream.P2.Temp;
        stream.P = SPLIT(stream.P2,stream.R2);
        [stream.FP,stream.F] = CalcFP(stream.F2,stream.P,-diff(QTEC));
        stream.R1 = SPLIT(stream.F2,stream.FP);        
        stmF0 = MIX(stream.R1,stream.F);
        TF0a = stmF0.Temp;
    end

    function dT = extTEHP1(x)
        dT = ones(size(x)); % 设初值
        stream.F0.Temp = x(1); stream.R2.Temp = x(2);
        stream.F0 = DCMD_PackStream(stream.F0);
        stream.R2 = DCMD_PackStream(stream.R2);
        [SOUTs,QTEC] = HXiTEHP([stream.F0,stream.R2],extTEC,opts,Eset);
        stream.F1 = SOUTs(1);
        stream.P1 = SOUTs(2);
        SINs = SOUTs;
        [profile,SOUTs] = TEHPiDCMD1(SINs,TECs,TEXs,membrane,flowPattern,opts);
        stream.F2 = SOUTs(1);
        stream.P2 = SOUTs(2);
        dT(2) = stream.P2.Temp-x(2);
        stream.P = SPLIT(stream.P2,stream.R2);
        [stream.FP,stream.F] = CalcFP(stream.F2,stream.P,-diff(QTEC));
        stream.R1 = SPLIT(stream.F2,stream.FP);        
        stream.F0 = MIX(stream.R1,stream.F);
        dT(1) = stream.F0.Temp-x(1);
    end

    function Eother = EDistributor(Ein,Eone)
        Eother = Ein-Eone;
    end

    function sout = Heater(sin,Q)
        sout = sin;
        sout.Enthalpy = sin.Enthalpy+Q;
        sout.Temp = sout.Enthalpy/sout.MassFlow/sout.SpecHeat;
    end

    function [stmFP,stmF,msg] = CalcFP(stmF2,stmP,E)
        stmFP = stmF2; stmF = stream.F; msg = ''; % 设初值
        % 联立HFP=E+HF-HP、HF=WF*cp*T0和HFP=(WF-WP)*cp*TF2求得WF
        stmF.MassFlow = (E-stmP.Enthalpy+stmP.MassFlow*stmF2.SpecHeat*stmF2.Temp)/stream.F.SpecHeat/(stmF2.Temp-stmF.Temp);
        stmF = DCMD_PackStream(stmF); % 维持进料温度为T0，计算stmF
        % 计算料液侧出料
        stmFP.MassFlow = stmF.MassFlow-stmP.MassFlow;
        stmFP.Temp = stmF2.Temp;
        stmFP = DCMD_PackStream(stmFP);
        dHB = E+stmF.Enthalpy-stmP.Enthalpy-stmFP.Enthalpy; % 系统能量平衡偏差
        if dHB > 1e-6
            msg = sprintf('【注意】CalcFP()能量平衡偏差%.5g！',dHB);
        end
    end

    function stmR2 = CalcR2(stmP1,stmP2)
        stmR2 = stmP2;
        stmR2.MassFlow = stmP1.MassFlow;
        stmR2 = DCMD_PackStream(stmR2);
    end
    
    function [sout1,sout2] = SPLIT(sin,var)
        sout1 = sin; % 初值
        if nargout == 2 && isa(var,'double') % 输入回流比计算2股出料
            R = var;
            sout1.MassFlow = sin.MassFlow/(R+1);
            sout2 = sin; % 设初值
            sout2.MassFlow = sin.MassFlow-sout1.MassFlow;
        else % 输入其中一股出料计算另一股出料
            sout2 = var;
            sout1.MassFlow = sin.MassFlow-sout2.MassFlow;
        end
        sout1 = DCMD_PackStream(sout1);
        sout2 = DCMD_PackStream(sout2);
    end

    function stmF = CalcFeed(stmP,stmFP)
        stmF = Stream;
        stmF.Temp = 298.15; % 进料温度固定为298.15K
        stmF.MassFlow = stmP.MassFlow+stmFP.MassFlow; % 进料为纯水
        stmF = DCMD_PackStream(stmF);
    end

    function sout = MIX(sin1,sin2)
        sout = sin1; % 设初值
        sout.MassFlow = sin1.MassFlow+sin2.MassFlow; % 质量守恒
        sout.Enthalpy = sin1.Enthalpy+sin2.Enthalpy; % 能量守恒
        sout.Temp = sout.Enthalpy/sout.MassFlow/sout.SpecHeat; 
    end

    function result = IsConvergence(sins1,sins2,tol)
        if ~exist('tol','var')
            tol = 1e-5;
        end
        if isstruct(sins1)
            vec1a = arrayfun(@(x)x.MassFlow,sins1);
            vec1a = vec1a/mean(vec1a);
            vec1b = arrayfun(@(x)x.Temp,sins1);
            vec1b = vec1b/mean(vec1b);
            vec1 = [diff(vec1a),diff(vec1b)];
        else
            vec1 = diff(sins1);
        end
        if isstruct(sins2)
            vec2a = arrayfun(@(x)x.MassFlow,sins2);
            vec2a = vec2a/mean(vec2a);
            vec2b = arrayfun(@(x)x.Temp,sins2);
            vec2b = vec2b/mean(vec2b);
            vec2 = [diff(vec2a),diff(vec2b)];
        else
            vec2 = diff(sins2);
        end
        d = norm([vec1,vec2]); % disp(d)
        result = (d<tol);
    end

end