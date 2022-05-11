function Postprocess(results,opStr)
%% 后处理
switch opStr
    case 'plotBar'
        % 对比外置集成半导体热泵与传统DCMD的单位能耗
        SEC0 = results.SEC(strcmp(results.CFG,'classical'));
        SEC1 = results.SEC(strcmp(results.CFG,'extTEHP'));
        bar([SEC0,SEC1])
        xlabel('Case No.#')
        ylabel('SEC [kWh/kg]')
        legend({'classical','extTEHP'})
    case 'RSM'
        % 响应面分析
        factors = [results.W1,results.T1,results.W2,results.T2];
        fStr = ['W1';'T1';'W2';'T2'];
        responses = [results.SEC];
        rStr = 'SEC';
        alpha = 0.01; % 显著性水平
        rstool(factors,responses,'quadratic',alpha,fStr,rStr)
end