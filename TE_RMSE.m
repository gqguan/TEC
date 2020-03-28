function [ RMSE output ] = TE_RMSE( x, TEC, ExpData )
%% calculate the RMSE between predicted and experimental results
%  notes of I/O arguments
%  x       - (i double array) [NumRatio GeomFactor] of thermocouples
%  TEC     - (i structure) initial parameters of thermocouples,
%                          ref. to "TE_ImportExpData.m"
%  ExpData - (i table) experimental results, ref. to "TE_ImportExpData.m"
%  RMSE    - (o double scalar) RMSE(ExpData.QH-QH)
%  output  - (o structure) .TEC: output TEC parameters
%                          .results: list results of exp. vs sim.
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-09
%
%% function body
% initialize
QH  = zeros(size(ExpData.QH));
QC  = zeros(size(ExpData.QC));
% reset TEC according to the given x
switch length(x)
    case 1 % x = GeomFactor
        TEC.GeomFactor = x;
    case 2 % x = [NumRatio GeomFactor]
        TEC.NumRatio   = x(1);
        TEC.GeomFactor = x(2);
end
%
% 计算理论吸、放热量
NumExpData = height(ExpData);
for i = 1: NumExpData
    % TEC冷热侧温度
    TC = ExpData.TC(i)+273.15; TH = ExpData.TH(i)+273.15;
    % 计算电流上下边界
    IBound = TE_Current(TH, TH, TEC); % 冷侧温度与热侧相同
    IMax = max(IBound);
    IMin = min(IBound);
    % 判定实验测得电流是否在理论范围
    if (ExpData.I(i) > IMax || ExpData.I(i) < IMin)
        fprintf('Given current %5.3f is out of range [%5.3f %5.3f]!\n', ...
                ExpData.I(i), IMin, IMax);
        output = [];
        RMSE = [];
        return;
    end
    % 输入实验时TEC的电压和电流
    TEC.Voltage = ExpData.U(i);
    TEC.Current = ExpData.I(i);
    % 计算TEC吸放热量
    [Q, TEC] = TE_Heat(TH, TC, TEC);
    QH(i) = Q(1);
    QC(i) = Q(2);
%     COP(i) = QC(i)./(QH(i)-QC(i));
end
% 输出结果
output.TEC = TEC;
tabout = table(ExpData.QH, 'VariableNames', {'QH_EXP'});
tabout = [tabout,table(ExpData.QC, 'VariableNames', {'QC_EXP'})];
tabout = [tabout,table(QH, 'VariableNames', {'QH_SIM'})];
tabout = [tabout,table(QC, 'VariableNames', {'QC_SIM'})];
output.results = tabout;
% COP_exp = ExpData.QC./(ExpData.QH-ExpData.QC);
RMSE = MVA_diff(ExpData.QH, QH, 'RMSE')/mean(ExpData.QH) ...
      +MVA_diff(ExpData.QC, QC, 'RMSE')/mean(ExpData.QC);
% RMSE = MVA_diff(ExpData.QH, QH, 'RMSE');
%
end

