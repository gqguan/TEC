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

%% 清除临时变量
clearvars data raw stringVectors;