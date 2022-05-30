%% 模拟一维集成热泵DCMD膜组件中的传热和传质现象：传统和集成半导体热泵DCMD系统的单位能耗
% 对比不同DCMD系统（classical、extTEHP、feedTEHP和permTEHP）给定输入能量的最大产水率
%
% by Dr. Guan Guoqiang @ SCUT on 2022/05/21

%% 初始化
clear
% 添加\Common目录以便调用自定义公用函数
homePath = cd;
idxPath = strfind(homePath,'\');
addpath([homePath(1:idxPath(2)),'Common\'])

%% 模拟DCMD系统能耗情况
CFGLvls = {'classical'}; %  DCMD配置方案
ELvls = num2cell(5:1.25/2:20);
% 根据各因素水平生成全因素实验方案
dMat = fullfact([length(CFGLvls),length(ELvls)]);
nLvl = max(dMat);
CaseNo = GenerateCaseNo('ffd',dMat);
% 实验条件
% dTab = table('size',size(dMat)+[0,1],...
%              'VariableTypes',{'cell','cell','double'},...
%              'VariableNames',{'SN','CFG','E'}); 
nCase = size(dMat,1);
SN = cell(1,nCase);
CFG = cell(1,nCase);
E = cell(1,nCase);
for i = 1:nCase
    SN(i) = CaseNo(i);
    CFG(i) = CFGLvls(dMat(i,1));
    E(i) = ELvls(dMat(i,2));
end
varNames = {'SN' 'CFG' 'W1' 'T1' 'W2' 'T2' 'RR' 'WF' 'WP' 'QM' 'Q1' 'E1'...
            'Q2' 'E2' 'QTEC' 'NTEC' 'TF0' 'SEC' 'NOTE' 'profile'};
% varTypes = {'cell' 'cell' 'double' 'double' 'double' 'double' 'double' ...
%             'double' 'double' 'double' 'double' 'double' 'double'...
%             'double' 'double' 'double' 'double' 'double' 'cell' 'struct'};
% results = table('size',[size(dMat,1),length(varNames)],...
%              'VariableTypes',varTypes,...
%              'VariableNames',varNames); 
results = table;
hbar = parfor_progressbar(nCase,'Computing...');
parfor iCase = 1:nCase % parfor循环中的变量为临时变量，不能在parfor循环以外访问
    [maxWP,x,tbl] = MaxWP(E{iCase},CFG{iCase});
    results(iCase,:) = tbl;
    hbar.iterate(1)
end
results.Properties.VariableNames = varNames;
close(hbar);

%% 后处理
% 响应面分析
% Postprocess(results,'plotBar');
