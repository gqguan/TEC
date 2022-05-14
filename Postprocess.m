function Postprocess(results,opStr)
colorCode = {'#0072BD' '#D95319' '#77AC30' '#A2142F' '#EDB120' '#4DBEEE' '#7E2F8E'};
%% 后处理
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
    case 'RSM'
        % 响应面分析
        factors = [results.W1,results.T1,results.W2,results.T2];
        fStr = ['W1';'T1';'W2';'T2'];
        responses = [results.SEC];
        rStr = 'SEC';
        alpha = 0.01; % 显著性水平
        rstool(factors,responses,'quadratic',alpha,fStr,rStr)
end