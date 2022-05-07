function [RR,QM,WF,WP,TP1,TP2] = CalcReflux(profile,Q1)
    % 已知膜组件温度分布，计算DCMD系统在给定料液加热量和进、出料流率和温度下的回流比
    if nargin ~= 2
        error('缺少CalcReflux()输入参数！')
    end
    T0 = 298.15; % 进料温度为环境温度[K]
    QM = sum(profile.QM); % 跨膜传热量[W]
    WP = sum(arrayfun(@(x)x.MassFlow,profile.SM)); % 跨膜渗透量[kg/s]
    TM = mean(profile.T(3:4,:),2);
    iStart = strfind(profile.Remarks,'：');
    % 膜组件料液侧的进口流率
    W1 = profile.S1(1).MassFlow;
    % 膜组件料液侧进、出口温度
    TF1 = profile.S1(1).Temp;
    TF2 = profile.S1(end).Temp;
    % 膜组件渗透侧进出口温度
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
    % 计算回流比的推导见case_DeriveFormulus.m
    TMF = TM(1);
    RR = -(Q1*TF2 - QM*T0 + T0*TF1*W1*cp1 - TF1*TF2*W1*cp1 + T0*TF2*WP*cp1 - T0*TMF*WP*cp1)/(TF2*(Q1 - QM + T0*WP*cp1 - TMF*WP*cp1));
    WF = (W1+RR*WP)/(1+RR);
end