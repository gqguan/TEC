%% show the results in both prediction and experiment with figures
% 
% by Dr. Guan Guoqiang @ SCUT on 2019-08-10
%
% 检查工作空间中TEC和ExpData变量是否存在
wrkvars_exist = exist('TEC', 'var') & exist('ExpData', 'var');
if wrkvars_exist
    QH = zeros(size(ExpData.QH));
    QC = zeros(size(ExpData.QC));
    for i = 1:height(ExpData)
        Q = TE_Heat(ExpData.TH(i), ExpData.TC(i), ExpData.I(i), TEC);
        QH(i) = Q(1);
        QC(i) = Q(2);
    end
    % 计算COP
    COP = QC./(QH-QC);
    COP_exp = ExpData.QC./(ExpData.QH-ExpData.QC);
    % 输出结果
    plot(ExpData.QC, COP_exp, 'ro', QC, COP, 'b*');
    xlabel('QC, [W]'); ylabel('COP'); 
    legend('Exp', 'Sim', 'Location', 'bestoutside');
    fprintf('RMSE.QH(exp-sim)/exp = %5.3f\n', ...
            norm((ExpData.QH-QH)./ExpData.QH));
    fprintf('RMSE.QC(exp-sim)/exp = %5.3f\n', ...
            norm((ExpData.QC-QC)./ExpData.QC));
else
    fprintf('Abort as the TEC and ExpData are not existed!\n');
end
clear wrkvars_exist i Q;