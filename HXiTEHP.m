function [SOUTs,QTEC] = HXiTEHP(SINs,TEC,opts,Eset)
% 给定冷热侧进料，其中热侧为SINs(1)、冷侧为SINs(2)及集成半导体热泵的操作参数
% 计算冷热侧出料，其中热侧为SOUTs(1)、冷侧为SOUTs(2)
% 若指定Eset则按Eset计算TEC功耗及相应的出料
    % 初始化
    strIU = {'Current','Voltage'};
    strUnit = {'A','V'};
    % 获取DuctGeom（流道几何尺寸）、Stream（物料定义）、MembrProps（膜材料性质）
    [~,Stream,MembrProps] = InitStruct();
    flowPattern = 'countercurrent';
    % 调试用输入参数
    if nargin == 0
        T1 = 273.15+50; T2 = 273.15+45; % [K]
        W1 = 1.217e-2; W2 = 6.389e-3; % [kg/s]
        s1 = Stream;
        s1.Temp = T1;
        s1.MassFlow = W1;
        SINs(1) = DCMD_PackStream(s1); % 膜组件料液侧进料
        s2 = s1;
        s2.Temp = T2;
        s2.MassFlow = W2;
        SINs(2) = DCMD_PackStream(s2); % 膜组件渗透侧进料
        iTEC = 22;
        % set properties for all TECs
        load('TEC_Params.mat','TEC_Params') % 载入已有的TEC计算参数
        TEC = TEC_Params.TEC(iTEC);
        % TEC计算参数
        opts = [TEC_Params.opt1(iTEC),TEC_Params.opt2(iTEC)];
        % 设定TEC功耗
        Eset = 12;
    end
    
    if exist('Eset','var')
        z0 = 1.5;
        [z,dE,exitflag] = fzero(@(z)CalcEC(z)-Eset,z0);
    else
        CalcEC()
    end
    
    
    function EC = CalcEC(argin)
        if nargin == 1
            TEC.(strIU{opts(2)+1}) = argin;
        end
        TF0 = SINs(1).Temp;
        TP2 = SINs(2).Temp;
        % 假定TEC冷热侧水箱加热过程为全混流，求满足能量平衡的冷热水箱温度
        solOpt = optimoptions(@lsqnonlin,'Display','none',...
                                     'Algorithm','trust-region-reflective');
        x0 = [TF0,TP2];
        lb = [298.15,273.15+5];
        ub = [273.15+95,273.15+95];
        [x,resnorm,~,exitflag] = lsqnonlin(@HeatBalance,x0,lb,ub,solOpt);
        EC = QTEC(1)-QTEC(2);
        if resnorm > 1e-4
            warning('【注意】CalcEC()求解外置集成半导体热泵加热和冷却单元进料温度分别为%.5gK和%.5gK，其结果偏差%.5g（exitflag=%d）',...
                x(1),x(2),resnorm,exitflag)
            return
        end
    end

    function dH = HeatBalance(T)
        SOUTs = SINs;
        SOUTs(1).Temp = T(1);
        SOUTs(2).Temp = T(2);
        % 计算出料焓值
        SOUTs(1) = DCMD_PackStream(SOUTs(1)); 
        SOUTs(2) = DCMD_PackStream(SOUTs(2));
        TH = T(1); TC = T(2);
        QTEC = TE_Heat(TH,TC,TEC,opts(1),opts(2));
        dH(1) = abs(SINs(1).Enthalpy+QTEC(1)-SOUTs(1).Enthalpy);
        dH(2) = abs(SINs(2).Enthalpy-QTEC(2)-SOUTs(2).Enthalpy);
    end

end