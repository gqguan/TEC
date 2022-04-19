function [ RMSE, output ] = TE_RMSE( x, TEC, ExpData, opt )
%% calculate the RMSE between predicted and experimental results
%  notes of I/O arguments
%  x       - (i double array) [NumRatio GeomFactor] of thermocouples
%  TEC     - (i structure) initial parameters of thermocouples,
%                          ref. to "TE_ImportExpData.m"
%  ExpData - (i table) experimental results, ref. to "TE_ImportExpData.m"
%  opt     - (i integer scalar) optional running mode
%            = 0 (default): 按参考文献[1]计算吸放热量
%            = 1          : 按参考文献[2]计算吸放热量
%  RMSE    - (o double scalar) RMSE
%  output  - (o structure) .TEC: output TEC parameters
%                          .results: list results of exp. vs sim.
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-09
%
%  References
%  [1] Xuan X C, et al. Cryogenics, 2002, 42: 273-278.
%  [2] Huang B J, et al. International Journal of Refrigeration 2000, 23(3): 208-218.
%
%% function body
% initialize
% 运行模式参数设定
switch nargin
    case(3) % 函数调用时不输入opt
        opt = 0;
    case(4)
        if opt ~= 0 && opt ~= 1
            prompt = sprintf('Unknown specified running mode of %d for TE_RMSE()', opt);
            TE_log(prompt, 1);
            return
        end
end
%
QH  = zeros(size(ExpData.QH));
QC  = zeros(size(ExpData.QC));
switch opt
    case(0)
        % reset TEC according to the given x
        switch length(x)
            case(1) % x = GeomFactor
                TEC.GeomFactor = x;
            case(2) % x = [NumRatio GeomFactor]
                TEC.NumRatio   = x(1);
                TEC.GeomFactor = x(2);
            case(3)
                TEC.NumRatio   = x(1);
                TEC.GeomFactor = x(2);
                TEC.HTCoefficient = x(3);
        end
    case(1) % TEC参数为定值
%         if length(x) ~= 3
%             prompt = sprintf('Incorrect Length(x) = %d', length(x));
%             TE_log(prompt, 1);
%             return
%         end
        % reset TEC parameters
        TEC.Parameters = x;
end
%
% 计算理论吸、放热量
NumExpData = height(ExpData);
for i = 1: NumExpData
    % TEC冷热侧温度
    TC = ExpData.TC(i)+273.15; TH = ExpData.TH(i)+273.15;
%     % 计算电流上下边界
%     IBound = TE_Current(TH, TH, TEC, opt); % 冷侧温度与热侧相同
%     IMax = max(IBound);
%     IMin = min(IBound);
%     % 判定实验测得电流是否在理论范围
%     if (ExpData.I(i) > IMax || ExpData.I(i) < IMin)
%         prompt = sprintf('Given current %5.3f is out of range [%5.3f %5.3f]', ExpData.I(i), IMin, IMax);
%         TE_log(prompt, 1);
%         output = [];
%         RMSE = [];
%         return;
%     end
    % 输入实验时TEC的电压和电流
    TEC.Voltage = ExpData.U(i);
    TEC.Current = ExpData.I(i);
    % 计算TEC吸放热量
    [Q, TEC] = TE_Heat(TH, TC, TEC, opt, 1);
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
RMSE = sqrt(sum(((tabout.QH_EXP-tabout.QH_SIM)./tabout.QH_EXP).^2 + ...
                ((tabout.QC_EXP-tabout.QC_SIM)./tabout.QC_EXP).^2));
% COP_exp = ExpData.QC./(ExpData.QH-ExpData.QC);
% RMSE = MVA_diff(ExpData.QH, QH, 'RMSE')/mean(ExpData.QH) ...
%       +MVA_diff(ExpData.QC, QC, 'RMSE')/mean(ExpData.QC);
% RMSE = MVA_diff(ExpData.QH, QH, 'RMSE');
%
end