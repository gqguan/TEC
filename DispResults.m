function p = DispResults(profile,DuctGeom,membrane)
    persistent iLine
    if isempty(iLine)
        iLine = 0;
    end
    iLine = iLine+1;
    % 产水量
    totalWP = sum([profile.SM.MassFlow]); % [kg/s]
    JM = totalWP/membrane.Area*3600; % [kg/m2-h]
    fprintf('产水量为%.4e[kg/s]（渗透通量为%.3f[kg/m2-h]）\n', totalWP, JM)
    % TEC特性
    TECE = sum(cell2mat(cellfun(@(x)x(:,1)-x(:,2),profile.QTEC,'UniformOutput',false)),2);
    Q1 = sum(cell2mat(cellfun(@(x)x(:,1),profile.QTEC,'UniformOutput',false)),2);
    Q2 = sum(cell2mat(cellfun(@(x)x(:,2),profile.QTEC,'UniformOutput',false)),2);
    fprintf('DCMD膜组件料液侧加热TEC电功率为%.4f[W]：放热量%.4f[W]，从环境吸热量%.4f[W]\n', ...
        TECE(1), Q1(1), Q2(1))
    fprintf('DCMD膜组件渗透侧冷却TEC电功率为%.4f[W]：吸热量%.4f[W]，向环境放热量%.4f[W]\n', ...
        TECE(2), Q2(2), Q1(2))
    % 绘制温度侧型
    lineColor = {'#7E2F8E';'#A2142F';'#D95319';'#77AC30';'#4DBEEE';'#0072BD'};
    lineStyle = {'-';'--';':';'-.'};
    x = linspace(0,DuctGeom.Length,length(profile.QM));
    p = plot([x;x;x;x;x;x]',profile.T',lineStyle{iLine});
    for iPlot = 1:length(p)
        p(iPlot).Color = lineColor{iPlot};
    end
    if iLine <= 4
        hold on
    else
        hold off
        iLine = 0;
    end
end
