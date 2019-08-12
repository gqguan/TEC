%% optimize geometry factor of thermocouple for 2-stage serial TEC
% 
% by Dr. Guan Guoqiang @ SCUT on 2019-08-11
%
clear;
% import data from text file "ExpData.txt"
fprintf('Import data from "ExpData.txt" ... \n');
tic
TE_ImportExpData
ExpData.TH = ExpData.TH+273.15;
ExpData.TC = ExpData.TC+273.15;
% initialize TEC parameters
TEC = struct('NumTC', 127, 'NumRatio', 0, 'GeomFactor', 8e-4, ...
             'HTCoefficient', 1200, 'HTArea', 0.0016);
% 计时
toc
% 求使计算最接近实验结果的热电偶几何因素 x = [HTCoefficient,GF]
fprintf('Calculating norm(COP(QC)_exp-COP(QC)_sim) ... \n');
x0 =  [TEC.HTCoefficient,TEC.GeomFactor];
[x, RMSE, exitflag] = fminsearch(@(x)TE_RMSE(x, TEC, ExpData), x0);
TEC.HTCoefficient = x(1);
TEC.GeomFactor    = x(2);
% 计时
toc
% 清理变量
clear N0;
% 输出
TE_ShowDiff;