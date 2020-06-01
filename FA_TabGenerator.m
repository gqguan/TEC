%% �������ر�
%
% �������˵��
% FactorNames - (i, string array) ������������
% LevelSets   - (i, double array) �����ص�ˮƽֵ
% FA_Tabout   - (o, table) ���ر�
% flag        - (o, integer) �������н��״̬������1Ϊ������0Ϊ�쳣
% 
% by Dr. Guan Guoqiang @ SCUT @ 2020/5/31
%
function [FA_Tabout,flag] = FA_TabGenerator(FactorNames,LevelSets)
%% ����������
% ����������������ر����������ˮƽ��
[NumFactor,NumLevel] = size(LevelSets);
if NumFactor ~= length(FactorNames)
    flag = 0;
    prompt1 = sprintf('Unmatched sizes of input arguements');
    disp(prompt1)
    return
end

%% ����ˮƽȫ������ƹ����ʼ���ر�
% �������ؾ���
Factor = ff2n(NumFactor);
% ���ɱ�
FA_Tabout = array2table(Factor,'VariableNames',FactorNames);

%% ������������ˮƽֵ��������ر�
for i = 1:NumFactor
    FA_Tabout{(FA_Tabout{:,i} == 0),i} = LevelSets(i,1);
    FA_Tabout{(FA_Tabout{:,i} == 1),i} = LevelSets(i,2);
end

%% �������״̬
flag = 1;
