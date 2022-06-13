%% 模拟给定功耗下冷热侧流率变化对集成半导体热泵DCMD系统能耗影响
% 模拟方案说明：
% 设定功耗为3水平，模拟冷热侧流率分别从Re=50~1000时集成半导体热泵DCMD系统能耗
%
% by Dr. Guan Guoqiang @ SCUT on 2022/06/13

%% 初始化
clear
% 添加\Common目录以便调用自定义公用函数
homePath = cd;
idxPath = strfind(homePath,'\');
addpath([homePath(1:idxPath(2)),'Common\'])

%% 模拟DCMD系统能耗情况
% 全因素设计
dMat = fullfact([1 3 3 3]);
nLvl = max(dMat);
CaseNo = GenerateCaseNo('ffd',dMat);
CFGLvls = {'extTEHP'}; %  DCMD配置方案
ELvls = {10,25,40}; % 设定功耗 [W]
W1Lvls = linspace(1.217e-4*5,1.217e-2,nLvl(3)); %  料液侧膜组件进料流率 [kg/s] Re=10*5~1000
W2Lvls = linspace(1.217e-4*5,1.217e-2,nLvl(4)); % 渗透侧膜组件进料流率 [kg/s] Re=10*5~1000
dTab = table; 
for i = 1:size(dMat,1)
    tt = table;
    tt.SN = CaseNo(i);
    tt.CFG = CFGLvls(dMat(i,1));
    tt.E = ELvls(dMat(i,2));
    tt.W1 = W1Lvls(dMat(i,3));
    tt.W2 = W2Lvls(dMat(i,4));
    dTab = [dTab;tt];
end
n = size(dMat,1);
results = table;
prompt = sprintf('开始%d个案例的计算',n);
hbar = parfor_progressbar(n,prompt);
parfor iExp = 1:n % parfor循环中的变量为临时变量，不能在parfor循环以外访问
    SN = CaseNo(iExp);
    CFG = CFGLvls(dMat(iExp,1));
    E = ELvls{dMat(iExp,2)};
    W1 = W1Lvls(dMat(iExp,3));
    W2 = W2Lvls(dMat(iExp,4));
    tab1 = [cell2table(SN),cell2table(CFG),table(E,W1,W2)];
    [stream,QTEC,profile,exitflag] = SimDCMD3(E,W1,W2,298.15);
    WP = stream.P.MassFlow;
    SEC = E/WP;
    results(iExp,:) = [tab1,table(WP,SEC,stream),table(QTEC),cell2table({profile})];
    hbar.iterate(1)
end
results.Properties.VariableNames = {'SN' 'CFG' 'E' 'W1' 'W2' 'WP' 'SEC' 'Stream' 'QTEC' 'profile'};
close(hbar);