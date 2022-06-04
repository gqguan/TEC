%% 模拟一维集成热泵DCMD膜组件中的传热和传质现象：传统和集成半导体热泵DCMD系统的单位能耗
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
% 添加\Common目录以便调用自定义公用函数
homePath = cd;
idxPath = strfind(homePath,'\');
addpath([homePath(1:idxPath(2)),'Common\'])

%% 模拟DCMD系统能耗情况
% 全因素设计
dMat = fullfact([4 3 3 3 3]);
nLvl = max(dMat);
CaseNo = GenerateCaseNo('ffd',dMat);
CFGLvls = {'classical','extTEHP','feedTEHP','permTEHP'}; %  DCMD配置方案
W1Lvls = linspace(1.217e-4*5,1.217e-2,nLvl(2)); %  料液侧膜组件进料流率 [kg/s] Re=10*5~1000
T1Lvls = linspace(273.15+50,273.15+65,nLvl(3)); % 料液侧膜组件进料温度 [K]
W2Lvls = linspace(1.217e-4*5,1.217e-2,nLvl(4)); % 渗透侧膜组件进料流率 [kg/s] Re=10*5~1000
T2Lvls = linspace(273.15+30,273.15+45,nLvl(5)); % 渗透侧膜组件进料温度 [K]
dTab = table; 
for i = 1:size(dMat,1)
    tt = table;
    tt.SN = CaseNo(i);
    tt.CFG = CFGLvls(dMat(i,1));
    tt.W1 = W1Lvls(dMat(i,2));
    tt.T1 = T1Lvls(dMat(i,3));
    tt.W2 = W2Lvls(dMat(i,4));
    tt.T2 = T2Lvls(dMat(i,5));
    dTab = [dTab;tt];
end
% 是否从存盘数据中继续计算
switch input('是[1]/否[0]是否从存盘数据中继续计算：')
    case 1
        load('dat_case9a.mat','results')
        fprintf('载入结果文件%d条记录\n',height(results))
        % 找出未完成计算的案例
        idx = cellfun(@(x)isempty(x)|isequal(x,'reset'),results.NOTE);
        plan = dTab(idx,:);
        n = height(plan);
        prompt = sprintf('继续剩余%d个案例的计算',n);
        hbar = parfor_progressbar(n,prompt);
        parfor iExp = 1:n
            SN = plan.SN(iExp);
            CFG = plan.CFG(iExp);
            W1 = plan.W1(iExp);
            T1 = plan.T1(iExp);
            W2 = plan.W2(iExp);
            T2 = plan.T2(iExp);
            tab1 = [cell2table(SN),cell2table(CFG),table(W1,T1,W2,T2)];
            [tab2,profile] = SimDCMD2(SN{1},W1,T1,W2,T2,CFG{1});
            tbl(iExp,:) = [tab1,tab2,cell2table({profile})];
            hbar.iterate(1)
        end
        tbl.Properties.VariableNames = {'SN' 'CFG' 'W1' 'T1' 'W2' 'T2' ...
           'RR' 'WF' 'WP' 'QM' 'Q1' 'E1' 'Q2' 'E2' 'QTEC' 'NTEC' 'TF0' 'SEC' 'NOTE' 'profile'};
        results(idx,:) = tbl;
        close(hbar);
    case 0 
        results = table;
        % 实验条件
        RR = inf; % 回流比
        n = size(dMat,1);
        prompt = sprintf('开始%d个案例的计算',n);
        hbar = parfor_progressbar(n,prompt);
        parfor iExp = 1:n % parfor循环中的变量为临时变量，不能在parfor循环以外访问
            SN = CaseNo(iExp);
            CFG = CFGLvls(dMat(iExp,1));
            W1 = W1Lvls(dMat(iExp,2));
            T1 = T1Lvls(dMat(iExp,3));
            W2 = W2Lvls(dMat(iExp,4));
            T2 = T2Lvls(dMat(iExp,5));
            tab1 = [cell2table(SN),cell2table(CFG),table(W1,T1,W2,T2)];
            [tab2,profile] = SimDCMD2(SN{1},W1,T1,W2,T2,CFG{1});
            results(iExp,:) = [tab1,tab2,cell2table({profile})];
            hbar.iterate(1)
        end
        results.Properties.VariableNames = {'SN' 'CFG' 'W1' 'T1' 'W2' 'T2' ...
            'RR' 'WF' 'WP' 'QM' 'Q1' 'E1' 'Q2' 'E2' 'QTEC' 'NTEC' 'TF0' 'SEC' 'NOTE' 'profile'};
        close(hbar);
    otherwise
        return
end

%% 后处理
% 响应面分析
Postprocess(results,'plotBar');
