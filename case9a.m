%% 模拟一维集成热泵DCMD膜组件中的传热和传质现象：传统DCMD系统的单位能耗
% 当膜组件集成半导体热泵TEHP时，可载入相应的TEC参数，例如用TEC_Params.mat中的H28、H05等
% 当未集成TEHP时，可载入TEC_Params.mat中的H00（近似绝热的边界条件：无电功输入且导热系数很小）
%
% by Dr. Guan Guoqiang @ SCUT on 2022/05/02
%
% 在稳态操作条件下，料液加热所需热量Q1为使膜组件料液侧出口温度升高到膜组件进口温度；
% 渗透液冷却所需冷量Q2为膜组件渗透侧出口温度降低到膜组件进口温度。
% 由此，料液加热器能耗按电功100%转化为热量计算，即为E1 = Q1；
% 渗透液冷却器能耗按TEC制冷耗电量计算，即为E2 = Q2/COP2
%
% 对于最低单位能耗的全回流稳态操作，料液加热量Q1 = QM+WF*CP*(TMH-T0)，详细推导见笔记2022/5/3；
% 渗透液冷却吸热量Q2 = QM+WP*CP*(TMC-TP2)
%
% 为确定DCMD系统吸放热量Q1和Q2，需要获得跨膜传热量和渗透量（QM和WP）的基础上进一步确定膜面温度（TMH和TMC）
% 方法1：分别列出WP和QM的计算式联立求解得TMH和TMC；
% 方法2：对模拟结果中的温度分布（即T(3)和T(4)）求平均值

%% 初始化
clear

%% 模拟DCMD系统能耗情况
% 采用4因素Box-Behnken设计
dBB = bbdesign(4)+2; % 2是为了将三个水平代号[-1 0 1]转换为索引
% frLvls = linspace(1.217e-5,2.434e-2,3); %  两侧膜组件进料流率 [kg/s] Re = 1~2000
frLvls = linspace(1.217e-4,1.217e-2,3); %  两侧膜组件进料流率 [kg/s] Re=10~1000
T1Lvls = linspace(273.15+45,273.15+60,3); % 料液侧膜组件进料温度 [K]
T2Lvls = linspace(273.15+5,273.15+20,3); % 渗透侧膜组件进料温度 [K]
results = table;
% 实验条件
n = size(dBB,1);
hbar = parfor_progressbar(n,'Computing...');
parfor iExp = 1:n
    W1 = frLvls(dBB(iExp,1));
    T1 = T1Lvls(dBB(iExp,2));
    W2 = frLvls(dBB(iExp,3));
    T2 = T2Lvls(dBB(iExp,2));
    tab1 = table(W1,T1,W2,T2);
    [tab2,profile] = SimDCMD(W1,T1,W2,T2);
    results(iExp,:) = [tab1,tab2,cell2table({profile})];
    hbar.iterate(1)
end
results.Properties.VariableNames = {'W1' 'T1' 'W2' 'T2' 'WP' 'QM' 'Q1' ...
    'E1' 'Q2' 'QTEC' 'E2' 'SEC' 'profile'};
close(hbar);
