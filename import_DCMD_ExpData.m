%% 从指定位置导入DCMD实验结果
%
% by Dr. Guan Guoqiang @ SCUT on 2020-05-24

%% 检查工作区中是否存在实验数据变量ExpData
if ~exist('ExpData', 'var') 
    fprintf('Existed DCMD experimental data are used.\n')
else
    return
end   

%% 从文件对话框选取导入的数据文件
[File, PathName] = uigetfile('*.*', '选取DCMD实验结果ExpData.xlsx ...', 'Multiselect', 'off');

%% 导入数据
filename = [PathName,File];
[~, ~, raw] = xlsread(filename,'summary','A3:X26');
stringVectors = string(raw(:,[1,8,23,24]));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,[2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22]);

%% 创建输出变量
data = reshape([raw{:}],size(raw));

%% 创建表
ExpData_raw = table;

%% 将导入的数组分配给列变量名称
ExpData_raw.Date = stringVectors(:,1);
ExpData_raw.STTime = data(:,1);
ExpData_raw.Heater_V = data(:,2);
ExpData_raw.Heater_I = data(:,3);
ExpData_raw.Fan_V = data(:,4);
ExpData_raw.Fan_I = data(:,5);
ExpData_raw.Pumps = data(:,6);
ExpData_raw.VarName8 = stringVectors(:,2);
ExpData_raw.QF_IN = data(:,7);
ExpData_raw.QP_IN = data(:,8);
ExpData_raw.VF_IN = data(:,9);
ExpData_raw.VP_IN = data(:,10);
ExpData_raw.TF_IN = data(:,11);
ExpData_raw.VarName14 = data(:,12);
ExpData_raw.TP_IN = data(:,13);
ExpData_raw.VarName16 = data(:,14);
ExpData_raw.TF_OUT = data(:,15);
ExpData_raw.VarName18 = data(:,16);
ExpData_raw.TP_OUT = data(:,17);
ExpData_raw.VarName20 = data(:,18);
ExpData_raw.WP = data(:,19);
ExpData_raw.VarName22 = data(:,20);
ExpData_raw.wF_IN = categorical(stringVectors(:,3));
ExpData_raw.Membrane = categorical(stringVectors(:,4));

%% 整理原始实验数据
ExpData = ExpData_raw(1:11,[9 10 11 12 13 15 17 19 21 23]);
% 将实验数据中的温度单位转换为K
ExpData.TF_IN = ExpData.TF_IN+273.15;
ExpData.TP_IN = ExpData.TP_IN+273.15;
ExpData.TF_OUT = ExpData.TF_OUT+273.15;
ExpData.TP_OUT = ExpData.TP_OUT+273.15;
% 将实验数据中的产水率单位转换为kg/s
ExpData.WP = ExpData.WP/1000;
    
%% 清除临时变量
clearvars ExpData_raw data raw stringVectors;