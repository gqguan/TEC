%% 生成因素表
%
% 输入输出说明
% FactorNames - (i, string array) 因数名称向量
% LevelSets   - (i, double array) 各因素的水平值
% FA_Tabout   - (o, table) 因素表
% flag        - (o, integer) 函数运行结果状态，其中1为正常，0为异常
% 
% by Dr. Guan Guoqiang @ SCUT @ 2020/5/31
%
function [FA_Tabout,flag] = FA_TabGenerator(FactorNames,LevelSets)
%% 输入参数检查
% 根据输入参数得因素表的因素数和水平数
[NumFactor,NumLevel] = size(LevelSets);
if NumFactor ~= length(FactorNames)
    flag = 0;
    prompt1 = sprintf('Unmatched sizes of input arguements');
    disp(prompt1)
    return
end

%% 按二水平全因素设计构造初始因素表
% 生成因素矩阵
Factor = ff2n(NumFactor);
% 生成表
FA_Tabout = array2table(Factor,'VariableNames',FactorNames);

%% 按各因素输入水平值代入的因素表
for i = 1:NumFactor
    FA_Tabout{(FA_Tabout{:,i} == 0),i} = LevelSets(i,1);
    FA_Tabout{(FA_Tabout{:,i} == 1),i} = LevelSets(i,2);
end

%% 输出运行状态
flag = 1;
