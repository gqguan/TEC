function [stream,QTEC,profile,exitflag] = SimDCMD3(Eset)
% 模拟给定功耗的外置集成半导体热泵DCMD系统
    % 初始化
    strIU = {'Current','Voltage'};
    strUnit = {'A','V'};
    % 获取DuctGeom（流道几何尺寸）、Stream（物料定义）、MembrProps（膜材料性质）
    [~,Stream,MembrProps] = InitStruct();
    flowPattern = 'countercurrent';
    % 设定外置集成半导体热泵换热单元的进料[stmF0,stmP2]初值
    if nargin == 0
        T1 = 273.15+50; T2 = 273.15+45; % [K]
        W1 = 1.217e-2; W2 = 6.389e-3; % [kg/s]
        s1 = Stream;
        s1.Temp = T1;
        s1.MassFlow = W1;
        stream.F0 = DCMD_PackStream(s1); % 膜组件料液侧进料
        s2 = s1;
        s2.Temp = T2;
        s2.MassFlow = W2;
        stream.R2 = DCMD_PackStream(s2); % 膜组件渗透侧进料
        iTEC1 = 22;
        % set properties for all TECs
        load('TEC_Params.mat','TEC_Params') % 载入已有的TEC计算参数
        extTEC = TEC_Params.TEC(iTEC1);
        opts = [TEC_Params.opt1(iTEC1),TEC_Params.opt2(iTEC1)];
        % 膜组件设定
        TECs(1) = TEC_Params.TEC(21); % 绝热边界
        TECs(2) = TEC_Params.TEC(21); % 绝热边界
        TEXs(1) = 298.15; % 膜组件的环境温度
        TEXs(2) = 298.15;
        membrane = MembrProps; % 膜性质
        % 设定系统操作参数
        Eset = 12; % TEC功耗
        refluxRatio = inf; % 回流比
    end
    % 计算外置集成半导体热泵的热、冷两侧进口流率和温度
%     x0 = [W1,T1+15,W2,T2+5];
%     solOpt = optimoptions(@lsqnonlin,'Display','iter',...
%                                  'Algorithm','levenberg-marquardt');
%     x1 = lsqnonlin(@extTEHP,x0,[],[],solOpt);
    % 迭代求解[stmF0,stmP2]
    for iter = 1:1000 % 最大迭代次数为1000
        [stmF0,stmR2] = extTEHP(stream.F0,stream.R2);
        if IsConvergence([stmF0,stmR2],[stream.F0,stream.R2])
            break
        else
            stream.F0 = stmF0;
            stream.R2 = stmR2;
        end
    end
    if iter == 1000
        exitflag = 0;
    else
        exitflag = 1;
    end
    outTab = table;
    stmList = {"F","FP","P","F0","F1","F2","R1","P1","P2","R2"};
    for iStm = 1:length(stmList)
        Stream = stmList{iStm};
        tab1r = [table(Stream),struct2table(stream.(stmList{iStm}))];
        outTab = [outTab;tab1r];
    end
    disp(outTab)
    
    function [stmF0,stmR2] = extTEHP(F0,R2)
        stream.F0 = F0; stream.R2 = R2;
        [SOUTs,QTEC] = HXiTEHP([stream.F0,stream.R2],extTEC,opts,Eset);
        stream.F1 = SOUTs(1);
        stream.P1 = SOUTs(2);
        SINs = SOUTs;
        [profile,SOUTs] = TEHPiDCMD1(SINs,TECs,TEXs,membrane,flowPattern,opts);
        stream.F2 = SOUTs(1);
        stream.P2 = SOUTs(2);
        stmR2 = CalcR2(stream.P1,stream.P2);
        stream.P = SPLIT(stream.P2,stream.R2);
        [stream.FP,stream.R1] = SPLIT(stream.F2,refluxRatio);
        stream.F = CalcFeed(stream.P,stream.FP);
        stmF0 = MIX(stream.R1,stream.F);
%         y(1) = (stmF0.MassFlow-x(1))/x(1); 
%         y(2) = (stmF0.Temp-x(2))/x(2);
%         y(3) = (stmR2.MassFlow-x(3))/x(3); 
%         y(4) = (stmR2.Temp-x(4))/x(4);
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
            tol = 1e-6;
        end
        vec1a = arrayfun(@(x)x.MassFlow,sins1);
        vec1a = vec1a/mean(vec1a);
        vec1b = arrayfun(@(x)x.Temp,sins1);
        vec1b = vec1b/mean(vec1b);
        vec1 = [vec1a,vec1b];
        vec2a = arrayfun(@(x)x.MassFlow,sins2);
        vec2a = vec2a/mean(vec2a);
        vec2b = arrayfun(@(x)x.Temp,sins2);
        vec2b = vec2b/mean(vec2b);
        vec2 = [vec2a,vec2b];
        d = norm(vec1-vec2); disp(d)
        result = (d<tol);
    end

end