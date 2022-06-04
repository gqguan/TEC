%% ģ��һά�����ȱ�DCMDĤ����еĴ��Ⱥʹ������󣺴�ͳ�ͼ��ɰ뵼���ȱ�DCMDϵͳ�ĵ�λ�ܺ�
% ��Ĥ������ɰ뵼���ȱ�TEHPʱ����������Ӧ��TEC������������TEC_Params.mat�е�H28��H05��
% ��δ����TEHPʱ��������TEC_Params.mat�е�H00�����ƾ��ȵı߽��������޵繦�����ҵ���ϵ����С��
%
% by Dr. Guan Guoqiang @ SCUT on 2022/05/02
%
% ����̬���������£���Һ������������Q1ΪʹĤ�����Һ������¶����ߵ�Ĥ��������¶ȣ�
% ��͸Һ��ȴ��������Q2ΪĤ�����͸������¶Ƚ��͵�Ĥ��������¶ȡ�
% �ɴˣ���Һ�������ܺİ��繦100%ת��Ϊ�������㣬��ΪE1 = Q1��
% ��͸Һ��ȴ���ܺİ�TEC����ĵ������㣬��ΪE2 = Q2/COP2
%
% ������͵�λ�ܺĵ�ȫ������̬��������Һ������Q1 = QM+WF*CP*(TMH-T0)����ϸ�Ƶ����ʼ�2022/5/3��
% ��͸Һ��ȴ������Q2 = QM+WP*CP*(TMC-TP2)
%
% Ϊȷ��DCMDϵͳ��������Q1��Q2����Ҫ��ÿ�Ĥ����������͸����QM��WP���Ļ����Ͻ�һ��ȷ��Ĥ���¶ȣ�TMH��TMC��
% ����1���ֱ��г�WP��QM�ļ���ʽ��������TMH��TMC��
% ����2����ģ�����е��¶ȷֲ�����T(3)��T(4)����ƽ��ֵ

%% ��ʼ��
clear
% ���\CommonĿ¼�Ա�����Զ��幫�ú���
homePath = cd;
idxPath = strfind(homePath,'\');
addpath([homePath(1:idxPath(2)),'Common\'])

%% ģ��DCMDϵͳ�ܺ����
% ȫ�������
dMat = fullfact([4 3 3 3 3]);
nLvl = max(dMat);
CaseNo = GenerateCaseNo('ffd',dMat);
CFGLvls = {'classical','extTEHP','feedTEHP','permTEHP'}; %  DCMD���÷���
W1Lvls = linspace(1.217e-4*5,1.217e-2,nLvl(2)); %  ��Һ��Ĥ����������� [kg/s] Re=10*5~1000
T1Lvls = linspace(273.15+50,273.15+65,nLvl(3)); % ��Һ��Ĥ��������¶� [K]
W2Lvls = linspace(1.217e-4*5,1.217e-2,nLvl(4)); % ��͸��Ĥ����������� [kg/s] Re=10*5~1000
T2Lvls = linspace(273.15+30,273.15+45,nLvl(5)); % ��͸��Ĥ��������¶� [K]
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
% �Ƿ�Ӵ��������м�������
switch input('��[1]/��[0]�Ƿ�Ӵ��������м������㣺')
    case 1
        load('dat_case9a.mat','results')
        fprintf('�������ļ�%d����¼\n',height(results))
        % �ҳ�δ��ɼ���İ���
        idx = cellfun(@(x)isempty(x)|isequal(x,'reset'),results.NOTE);
        plan = dTab(idx,:);
        n = height(plan);
        prompt = sprintf('����ʣ��%d�������ļ���',n);
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
        % ʵ������
        RR = inf; % ������
        n = size(dMat,1);
        prompt = sprintf('��ʼ%d�������ļ���',n);
        hbar = parfor_progressbar(n,prompt);
        parfor iExp = 1:n % parforѭ���еı���Ϊ��ʱ������������parforѭ���������
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

%% ����
% ��Ӧ�����
Postprocess(results,'plotBar');
