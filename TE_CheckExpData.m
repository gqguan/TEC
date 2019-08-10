%% check required experimental results existed in workspace
% 
% by Dr. Guan Guoqiang @ SCUT on 2019-08-10
%
% check ExpData existed or not
ExpData_exist = exist('ExpData', 'var');
if ExpData_exist
    fprintf('"ExpData" is existed. \n');
else
    % import data from text file "ExpData.txt"
    fprintf('Import data from "ExpData.txt" ... \n');
    TE_ImportExpData
    ExpData.TH = ExpData.TH+273.15;
    ExpData.TC = ExpData.TC+273.15;
end
clear ExpData_exist;