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
        Q = TE_Heat(ExpData.TH(i), ExpData.TC(i), ExpData.I(i), ...
                    TEC.NumTC, TEC.NumRatio, TEC.GeomFactor);
        QH(i) = Q(1);
        QC(i) = Q(2);
    end
    COP = QC./(QH-QC);
    COP_sim = ExpData.QC./(ExpData.QH-ExpData.QC);
    plot(QC, COP, 'ro', ExpData.QC, COP_sim, '*');
else
    fprintf('Abort as the TEC and ExpData are not existed!\n');
end
clear wrkvars_exist i Q;