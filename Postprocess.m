function results = Postprocess(results,opStr)
colorCode = {'#0072BD' '#D95319' '#77AC30' '#A2142F' '#EDB120' '#4DBEEE' '#7E2F8E'};
colorCode1 = {'#A2142F' '#D95319' '#EDB120' '#77AC30' '#4DBEEE' '#0072BD' ...
              '#A2142F' '#D95319' '#EDB120' '#77AC30' '#4DBEEE' '#0072BD' ...
              '#A2142F' '#D95319' '#EDB120' '#77AC30' '#4DBEEE' '#0072BD' ...
              '#A2142F' '#D95319' '#EDB120' '#77AC30' '#4DBEEE' '#0072BD'};
lineStyle1 = {'-' '-' '-' '-' '-' '-' ...
              '--' '--' '--' '--' '--' '--' ...
              '-.' '-.' '-.' '-.' '-.' '-.' ...
              ':' ':' ':' ':' ':' ':'};
%% 后处理
% 删除备注“注意”的数据
warnIdx = cellfun(@(x)contains(x,'【注意】'),results.NOTE);
switch input('需注意的存盘数据处理方式：[1]重置该记录以备重算；[2]删除相应的数据组：')
    case 1
        resetIdx = warnIdx;
        disp(results(resetIdx,:))
        results.NOTE(resetIdx) = {'reset'};
        fprintf('重置原数据集中%d个记录的NOTE为空值！\n',sum(resetIdx))
        return
    case 2
        grpId = arrayfun(@(x)ceil(x/4),1:height(results));
        idx = false(size(grpId));
        delGrpId = unique(grpId(warnIdx));
        for i = 1:length(delGrpId)
            idx = idx|(grpId==delGrpId(i));
        end
        disp(results(idx,:))
        results(idx,:) = [];
        fprintf('删除原数据集中%d组记录！\n',length(delGrpId))
    otherwise
end
%
switch opStr
    case 'plotBar'
        % 对比集成半导体热泵与传统DCMD的单位能耗
        cfgList = categories(categorical(results.CFG));
        % 绘制能耗比柱状图
        refCFG = 'classical';
        SEC0 = results.SEC(strcmp(results.CFG,refCFG));
        cfgList(strcmp(cfgList,refCFG)) = [];
        SEC = cell2mat(cellfun(@(x)results.SEC(strcmp(results.CFG,x)),cfgList','UniformOutput',false));
        RSEC = SEC./SEC0;
        b = bar(RSEC,'BaseValue',1);
        for i = 1:length(b)
            b(i).FaceColor = colorCode{i};
        end
        xlabel('Case No.#')
        ylabel('Relative SEC')
        legend(cfgList)
        fprintf('采用extTEHP、feedTEHP和permTEHP的平均能耗比分别为%.4g、%.4g和%.4g\n',mean(RSEC,'omitnan'))
        fprintf('最低能耗分别为%.4g、%.4g和%.4g\n',min(RSEC,[],'omitnan'))
    case 'TProfile' % 绘制膜组件中的温度侧形
        snList = {14 15 16};
        iLineMarked = ones(1,length(snList));
        skipLineNum = 6;
        for i = 1:length(snList), iLineMarked(i) = iLineMarked(i)+(i-1)*skipLineNum; end
        TProfile = cell2mat(cellfun(@(x)results.profile(x,1).T',snList,'UniformOutput',false));
        p = plot(TProfile);
        xlim([1,size(TProfile,1)-1])
        xticks([])
        for i=1:length(p), p(i).Color = colorCode1{i}; p(i).LineStyle = lineStyle1{i}; p(i).LineWidth = 1.5; end
        legendStr = cellfun(@(x)strcat(regexp(results.SN{x},'ffd\(\w*\)\|','match'),results.CFG(x)),snList);
        legend(p(iLineMarked),legendStr)
        xlabel('Length along DCMD module')
        ylabel('Temperature (K)')
        hold off
    case 'TEC' 
        snList = {13 14 15 16};
        for i = 1:length(snList)
            iSn = snList{i};
            cfg = results.CFG{iSn};
            switch cfg
                case 'classical'
                    ETEC = results.E2(iSn);
                    iStart = strfind(results.profile(iSn).Remarks,'：');
                    switch results.profile(iSn).Remarks(iStart+1:end)
                        case('cocurrent')
                            TP1 = results.profile(iSn).S2(1).Temp;
                            TP2 = results.profile(iSn).S2(end).Temp;
                        case('countercurrent')
                            TP1 = results.profile(iSn).S2(end).Temp;
                            TP2 = results.profile(iSn).S2(1).Temp;
                    end
                    TH = 298.15;
                    TC = mean([TP1 TP2]);
                    Q2TEC = results.Q2(iSn);
                    Q1TEC = Q2TEC+ETEC;
                case 'extTEHP'
                    ETEC = results.E2(iSn);
                    TF0 = results.TF0(iSn);
                    TF1 = results.profile(iSn).S1(1).Temp;
                    TF2 = results.profile(iSn).S2(end).Temp;
                    iStart = strfind(results.profile(iSn).Remarks,'：');
                    switch results.profile(iSn).Remarks(iStart+1:end)
                        case('cocurrent')
                            TP1 = results.profile(iSn).S2(1).Temp;
                            TP2 = results.profile(iSn).S2(end).Temp;
                        case('countercurrent')
                            TP1 = results.profile(iSn).S2(end).Temp;
                            TP2 = results.profile(iSn).S2(1).Temp;
                    end
                    TH = mean([TF0 TF1]);
                    TC = mean([TP1 TP2]);
                    Q1TEC = results.Q1(iSn);
                    Q2TEC = results.Q2(iSn);
                case 'feedTEHP'
                    ETEC = results.E2(iSn);
                    % TEC(1)热侧平均温度及放热量
                    TH = mean(results.profile(iSn).T(1,:));
                    Q1TEC = sum(cellfun(@(x)x(1,1),results.profile(iSn).QTEC));
                    % TEC(1)冷侧平均温度及吸热量
                    iStart = strfind(results.profile(iSn).Remarks,'：');
                    switch results.profile(iSn).Remarks(iStart+1:end)
                        case('cocurrent')
                            TP1 = results.profile(iSn).S2(1).Temp;
                            TP2 = results.profile(iSn).S2(end).Temp;
                        case('countercurrent')
                            TP1 = results.profile(iSn).S2(end).Temp;
                            TP2 = results.profile(iSn).S2(1).Temp;
                    end
                    TC = mean([TP1 TP2]);
                    Q2TEC = Q1TEC-ETEC;
                case 'permTEHP'
                    ETEC = results.E2(iSn);
                    % TEC(2)冷侧平均温度及吸热量
                    TC = mean(results.profile(iSn).T(6,:));
                    Q2TEC = sum(cellfun(@(x)x(2,2),results.profile(iSn).QTEC));
                    % TEC(2)热侧平均温度及放热量
                    TF0 = results.TF0(iSn);
                    TF1 = results.profile(iSn).S1(1).Temp;
                    TF2 = results.profile(iSn).S2(end).Temp;
                    iStart = strfind(results.profile(iSn).Remarks,'：');
                    switch results.profile(iSn).Remarks(iStart+1:end)
                        case('cocurrent')
                            TP1 = results.profile(iSn).S2(1).Temp;
                            TP2 = results.profile(iSn).S2(end).Temp;
                        case('countercurrent')
                            TP1 = results.profile(iSn).S2(end).Temp;
                            TP2 = results.profile(iSn).S2(1).Temp;
                    end
                    TH = mean([TF0 TF1]);
                    Q1TEC = Q2TEC+ETEC;
            end
            fprintf('%s：TEC热侧平均温度为%.4g[K]，放热量为%.4g[W]；',cfg,TH,Q1TEC)
            fprintf('冷侧平均温度为%.4g[K]，吸热量为%.4g[W]；',TC,Q2TEC)
            fprintf('平均温差%.4g[K]；COP2=%.4g\n',TH-TC,Q2TEC/ETEC)
        end
    case 'RSM'
        % 响应面分析
        factors = [results.W1,results.T1,results.W2,results.T2];
        fStr = ['W1';'T1';'W2';'T2'];
        responses = [results.SEC];
        rStr = 'SEC';
        alpha = 0.01; % 显著性水平
        rstool(factors,responses,'quadratic',alpha,fStr,rStr)
end