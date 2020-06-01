%% 生成因素表
%
% 输入输出说明
% SettingTab - (i, table) 各因数的水平设置表
% FA_Tabout  - (o, table) 因素表
% 
% by Dr. Guan Guoqiang @ SCUT @ 2020/6/1
%
function [FA_Tabout,flag] = FA_TabGenerator(SettingTab)
%% 检查输入参数
if size(SettingTab,1) ~= 1
    flag = 0;
    prompt1 = sprintf('[ERROR] Input table must have one row only!');
    disp(prompt1)
    return
end

%% 生成通用因素表
% 获取表头
FactorNames = SettingTab.Properties.VariableNames;
% 因素数
NumFactor = length(FactorNames);
% 获得各因素的水平设置数向量
LevelSettings = zeros(1, NumFactor);
for i = 1:NumFactor
    LevelSettings(i) = length(SettingTab{1,i});
end
% 生成因素矩阵
FactorMatrix = fullfact(LevelSettings);
% 生成表
FA_Tabout = array2table(FactorMatrix,'VariableNames',FactorNames);

%% 按各因素输入水平值代入的因素表
for i = 1:NumFactor
    LevelSets = SettingTab{1,i};
    for j = 1:length(LevelSets)
        FA_Tabout{(FA_Tabout{:,i} == j),i} = LevelSets(1,j);
    end
end

%% 输出运行状态
flag = 1;
