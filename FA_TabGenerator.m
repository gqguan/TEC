%% �������ر�
%
% �������˵��
% SettingTab - (i, table) ��������ˮƽ���ñ�
% FA_Tabout  - (o, table) ���ر�
% 
% by Dr. Guan Guoqiang @ SCUT @ 2020/6/1
%
function [FA_Tabout,flag] = FA_TabGenerator(SettingTab)
%% ����������
if size(SettingTab,1) ~= 1
    flag = 0;
    prompt1 = sprintf('[ERROR] Input table must have one row only!');
    disp(prompt1)
    return
end

%% ����ͨ�����ر�
% ��ȡ��ͷ
FactorNames = SettingTab.Properties.VariableNames;
% ������
NumFactor = length(FactorNames);
% ��ø����ص�ˮƽ����������
LevelSettings = zeros(1, NumFactor);
for i = 1:NumFactor
    LevelSettings(i) = length(SettingTab{1,i});
end
% �������ؾ���
FactorMatrix = fullfact(LevelSettings);
% ���ɱ�
FA_Tabout = array2table(FactorMatrix,'VariableNames',FactorNames);

%% ������������ˮƽֵ��������ر�
for i = 1:NumFactor
    LevelSets = SettingTab{1,i};
    for j = 1:length(LevelSets)
        FA_Tabout{(FA_Tabout{:,i} == j),i} = LevelSets(1,j);
    end
end

%% �������״̬
flag = 1;
