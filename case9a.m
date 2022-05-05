%% ģ��һά�����ȱ�DCMDĤ����еĴ��Ⱥʹ������󣺴�ͳDCMDϵͳ�ĵ�λ�ܺ�
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
% ����4����Box-Behnken���
dMat = bbdesign(4)+2; % 2��Ϊ�˽�����ˮƽ����[-1 0 1]ת��Ϊ����
% % ���ĸ������
% dMat = ccdesign(4)+3; % 2��Ϊ�˽�5��ˮƽ����[-2 -1 0 1 2]ת��Ϊ����
% % ȫ�������
% dMat = fullfact([3 3 3 3]);
nLvl = max(dMat);
% frLvls = linspace(1.217e-5,2.434e-2,3); %  ����Ĥ����������� [kg/s] Re = 1~2000
W1Lvls = linspace(1.217e-4,1.217e-2,nLvl(1)); %  ��Һ��Ĥ����������� [kg/s] Re=10~1000
T1Lvls = linspace(273.15+45,273.15+60,nLvl(2)); % ��Һ��Ĥ��������¶� [K]
W2Lvls = linspace(1.217e-4,1.217e-2,nLvl(1)); % ��͸��Ĥ����������� [kg/s] Re=10~1000
T2Lvls = linspace(273.15+5,273.15+20,nLvl(4)); % ��͸��Ĥ��������¶� [K]
results = table;
% ʵ������
n = size(dMat,1);
hbar = parfor_progressbar(n,'Computing...');
parfor iExp = 1:n % parforѭ���еı���Ϊ��ʱ������������parforѭ���������
    W1 = W1Lvls(dMat(iExp,1));
    T1 = T1Lvls(dMat(iExp,2));
    W2 = W2Lvls(dMat(iExp,3));
    T2 = T2Lvls(dMat(iExp,4));
    tab1 = table(W1,T1,W2,T2);
    [tab2,profile] = SimDCMD(W1,T1,W2,T2);
    results(iExp,:) = [tab1,tab2,cell2table({profile})];
    hbar.iterate(1)
end
results.Properties.VariableNames = {'W1' 'T1' 'W2' 'T2' 'WP' 'QM' 'Q1' ...
    'E1' 'Q2' 'QTEC' 'E2' 'NTEC' 'SEC' 'profile'};
close(hbar);

%% ����
% ��Ӧ�����
Postprocess(results);
