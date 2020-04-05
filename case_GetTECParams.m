%% 从半导体制冷片的性能测试实验结果中拟合TEC参数
%
% by Dr. GUAN Guoqiang @ SCUT on 2020-03-27
%
%% 初始化
clear;
% 数据结构定义
TEC = struct('NumTC', 190, 'NumRatio', 0.9, 'GeomFactor', 0.7e-3, ...
             'HTCoefficient', 270, 'HTArea', 0.0016, ...
             'SeebeckCoefficient', [], 'ElecConductance', [], ...
             'ThermConductance', [], 'Voltage', [], 'Current', [], ...
             'Parameters', []);
% 半导体制冷片的性能测试实验数据
% 从实验数据文件ExpData.txt中导入，实验数据存于工作空间的表变量ExpData中
TE_ImportExpData
%
%% 优化
% 参数初值
x0 = [TEC.NumRatio,TEC.GeomFactor];
x1 = [1,1,1;1,1,1;1,1,1];
% 定义优化参数
options = optimset('PlotFcns', @optimplotfval);
% 定义目标函数
f0 = @(x)(TE_RMSE(x, TEC, ExpData, 0));
f1 = @(x)(TE_RMSE(x, TEC, ExpData, 1));
% 获得优化参数
x = fminsearch(f1, x1, options);
%% 输出结果
[~,output] = f1(x);
% 输入TEC部件号
output.pid = input('Input TEC part no.: ', 's');
% 构造表
current_tab = struct2table(output, 'AsArray', 1);
% 当前目录存在TEC参数文件时载入表TEC_Params
if exist('TEC_Params.mat', 'file') == 2
    load('TEC_Params.mat')
    TEC_Params = [TEC_Params;current_tab];
else
    TEC_Params = current_tab;
end
% 结果存盘
save('TEC_Params.mat', 'TEC_Params')