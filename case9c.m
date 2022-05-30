%% ģ��һά�����ȱ�DCMDĤ����еĴ��Ⱥʹ������󣺴�ͳ�ͼ��ɰ뵼���ȱ�DCMDϵͳ�ĵ�λ�ܺ�
% �ԱȲ�ͬDCMDϵͳ��classical��extTEHP��feedTEHP��permTEHP��������������������ˮ��
%
% by Dr. Guan Guoqiang @ SCUT on 2022/05/21

%% ��ʼ��
clear
% ���\CommonĿ¼�Ա�����Զ��幫�ú���
homePath = cd;
idxPath = strfind(homePath,'\');
addpath([homePath(1:idxPath(2)),'Common\'])

%% ģ��DCMDϵͳ�ܺ����
CFGLvls = {'classical'}; %  DCMD���÷���
ELvls = num2cell(5:1.25/2:20);
% ���ݸ�����ˮƽ����ȫ����ʵ�鷽��
dMat = fullfact([length(CFGLvls),length(ELvls)]);
nLvl = max(dMat);
CaseNo = GenerateCaseNo('ffd',dMat);
% ʵ������
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
parfor iCase = 1:nCase % parforѭ���еı���Ϊ��ʱ������������parforѭ���������
    [maxWP,x,tbl] = MaxWP(E{iCase},CFG{iCase});
    results(iCase,:) = tbl;
    hbar.iterate(1)
end
results.Properties.VariableNames = varNames;
close(hbar);

%% ����
% ��Ӧ�����
% Postprocess(results,'plotBar');
