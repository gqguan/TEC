%% optimize geometry factor of thermocouple for 2-stage serial TEC
% 
% by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
clear;
% import data from text file "ExpData.txt"
fprintf('Import data from "ExpData.txt" ... \n');
tic
TE_ImportExpData
ExpData.TH = ExpData.TH+273.15;
ExpData.TC = ExpData.TC+273.15;
% initialize TEC parameters
TEC = struct('NumTC', 190, 'NumRatio', 1, 'GeomFactor', 3.8e-4);
% 计时
toc
% calculate the RMSE of QC
fprintf('Calculating RMSE of dQH(exp-sim) ... \n');
% 求使计算QC最接近实验结果的热电偶几何因素GF
[GF, RMSE, exitflag] = fminsearch(@(GF)TE_RMSE(GF, TEC, ExpData), 3.8e-4);
TEC.GeomFactor = GF;
% 计时
toc
% 清理变量
clear N0;
% 输出
fprintf('RMSE of Qc = %5.3f\n', RMSE);
TE_ShowDiff;