function Postprocess(results)
%% 后处理
% 响应面分析
factors = [results.W1,results.T1,results.W2,results.T2];
fStr = ['W1';'T1';'W2';'T2'];
responses = [results.SEC];
rStr = 'SEC';
alpha = 0.01; % 显著性水平
rstool(factors,responses,'quadratic',alpha,fStr,rStr)