function [Q,QM,WF,WP,TP1,TP2,TF1,TF2,dQ,TF0] = CalcHeat(profile,R,cfg)
    % 已知膜组件温度分布，计算DCMD系统在给定回流比时的料液加热量和进、出料流率
    if ~exist('R','var')
        R = inf; % 回流比为无穷大，即全回流
    end
    if ~exist('cfg','var')
        cfg = 'classical';
    end
    %
    T0 = 298.15; % 进料温度为环境温度[K]
    QM = sum(profile.QM); % 跨膜传热量[W]
    WP = sum(arrayfun(@(x)x.MassFlow,profile.SM)); % 跨膜渗透量[kg/s]
%     TM = sum((profile.T(3:4,1)+0.5*diff(profile.T(3:4,:),1,2)).*([1;1]*[profile.SM.MassFlow]),2)/WP;
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
            W2 = profile.S2(1).MassFlow;
        case('countercurrent')
            TP1 = profile.S2(end).Temp;
            TP2 = profile.S2(1).Temp;
            W2 = profile.S2(end).MassFlow;
        otherwise
            error('CalcHeat()输入参数profile字段Remarks中无有效的流型信息')
    end
    cp1 = mean([profile.S1.SpecHeat]);
    cp2 = mean([profile.S2.SpecHeat]);
    
    if isinf(R) % 全回流
        WF = WP;
        Q(1) = WF*cp1*(TF2-T0)+W1*cp1*(TF1-TF2); % 参见2022/5/16笔记，避免使用不精确的TM
        Q(2) = W2*cp2*(TP2-TP1); % 避免使用不精确的TM
        dQ = 0;
        TF0 = -(Q(1) - TF1*W1*cp1)/(W1*cp1); % 详见case_DeriveFormulus.m
    else % 部分回流，Q1的推导见case_DeriveFormulus.m
%         TMF = TM(1);
        WF = (W1+R*WP)/(1+R);
        switch cfg
            case {'classical','extTEHP'}
                % 料液侧加热量
                WR1 = R*(WF-WP);
                Q(1) = W1*cp1*TF1-(WF*cp1*T0+WR1*cp1*TF2); % 避免使用不精确的TM
%                 Q(1) = (QM*T0 + QM*R*TF2 - T0*TF1*W1*cp1 + TF1*TF2*W1*cp1 - T0*TF2*WP*cp1 + T0*TMF*WP*cp1 - R*T0*TF2*WP*cp1 + R*TF2*TMF*WP*cp1)/(TF2 + R*TF2);
                % 渗透液吸热量
                Q(2) = W2*cp2*(TP2-TP1); % 避免使用不精确的TM
%                 Q(2) = QM+WP*cp2*(TM(2)-TP2);
                dQ = 0;
                TF0 = -(Q(1) - TF1*W1*cp1)/(W1*cp1); % 详见case_DeriveFormulus.m
            case 'feedTEHP'
                % TEC传热量
                Q1 = sum(cellfun(@(x)x(1,1),profile.QTEC)); % 膜组件热侧TEC向料液放热量，若热侧集成TEHP放热量=加热料液的热量
                Q(1) = Q1;
                Q2 = sum(cellfun(@(x)x(1,2),profile.QTEC)); % 膜组件热侧TEC从环境吸热量，若热侧集成TEHP吸热量=冷却渗透液的冷量
                % 渗透液吸热量
                Q(2) = W2*cp2*(TP2-TP1); % 避免使用不精确的TM
%                 Q(2) = QM+WP*cp2*(TM(2)-TP2);
%                 Q(1) = (QM*T0 + QM*R*TF2 - T0*TF1*W1*cp1 + TF1*TF2*W1*cp1 - T0*TF2*WP*cp1 + T0*TMF*WP*cp1 - R*T0*TF2*WP*cp1 + R*TF2*TMF*WP*cp1)/(T0 + R*TF2);
%                 Q(1) = (W1-WP)*cp1*TF2+WP*cp1*TMF+QM-W1*cp1*TF1;
                WF = (TF1*W1 + R*TF2*WP)/(T0 + R*TF2);
                dQ2 = Q(2)-Q2;
                dQ = dQ2;
                TF0 = nan;
            case {'permTEHP','permTEHP1'}
                % TEC传热量
                Q1 = sum(cellfun(@(x)x(2,1),profile.QTEC)); % 膜组件热侧TEC向料液加热量
                Q2 = sum(cellfun(@(x)x(2,2),profile.QTEC));
                % 渗透液吸热量
                Q(2) = Q2;
%                 Q(2) = QM+WP*cp2*(TM(2)-TP2);
                WR1 = R*(WF-WP);
                Q(1) = W1*cp1*TF1-(WF*cp1*T0+WR1*cp1*TF2); % 避免使用不精确的TM
                TF0 = (WR1*cp1*TF2+WF*cp1*T0)/W1/cp1;
                % 详见case_DeriveFormulus.m
%                 Q(1) = (QM*T0 + QM*R*TF2 - T0*TF1*W1*cp1 + TF1*TF2*W1*cp1 - T0*TF2*WP*cp1 + T0*TMF*WP*cp1 - R*T0*TF2*WP*cp1 + R*TF2*TMF*WP*cp1)/(TF2 + R*TF2);
%                 TF0 = -(QM*T0 + QM*R*TF2 - T0*TF1*W1*cp1 - T0*TF2*WP*cp1 + T0*TMF*WP*cp1 - R*TF1*TF2*W1*cp1 - R*T0*TF2*WP*cp1 + R*TF2*TMF*WP*cp1)/(W1*cp1*(TF2 + R*TF2));
                dQ1 = Q(1)-Q1;
                dQ = dQ1;
            case 'permTEHP2'
                % TEC传热量
                Q1 = sum(cellfun(@(x)x(2,1),profile.QTEC)); % 膜组件热侧TEC向料液加热量
                Q2 = sum(cellfun(@(x)x(2,2),profile.QTEC));
                % 渗透液吸热量
                Q(2) = Q2;
%                 Q(2) = QM+WP*cp2*(TM(2)-TP2);
                WR1 = R*(WF-WP);
                Q(1) = W1*cp1*T1-WF*cp1*T0-WR1*cp1*TF2;
                TF0 = (W1*cp1*T1-WF*cp1*T0)/WR1/cp1;
                % 详见case_DeriveFormulus.m
%                 Q(1) = (QM*T0 + QM*R*TF2 - T0*TF1*W1*cp1 + TF1*TF2*W1*cp1 - T0*TF2*WP*cp1 + T0*TMF*WP*cp1 - R*T0*TF2*WP*cp1 + R*TF2*TMF*WP*cp1)/(TF2 + R*TF2);
%                 TF0 = -(QM*T0 - T0*TF1*W1*cp1 + TF1*TF2*W1*cp1 - T0*TF2*WP*cp1 + T0*TMF*WP*cp1 + R*TF1*TF2*W1*cp1 - R*T0*TF2*WP*cp1)/(R*(QM - TF1*W1*cp1 + TMF*WP*cp1));
                dQ1 = Q(1)-Q1;
                dQ = dQ1;
        end
        
    end
    
end