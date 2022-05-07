function maxQ2 = TE_ShowPerformance(Th,Tc,TEC,TECopts,silenceOrNot)
    % 日志文件
    logFilename = [cd,'\',datestr(datetime,'yyyymmdd'),'.log'];
    logger = log4m.getLogger(logFilename);
    logger.setFilename(logFilename);
    % 初始化
    if ~exist('silenceOrNot','var')
        silenceOrNot = true;
    end
    % 函数调试用
    if nargin == 0
        Th = 298.15; % [K]
        Tc = 278.15; % [K]
        load('TEC_Params.mat', 'TEC_Params')
        TEC = TEC_Params.TEC(14,:);
        TECopts = [1,0];
        silenceOrNot = false;
    end
    % 设定电流范围
    lb = 0.5; ub = 10; % 电流范围0.5~10[A]
    % 计算最大吸热量
    fun = @(x)-GetTECHeat(x,'cooling',Th,Tc,TEC,TECopts); % 负值求最小即求最大吸热量
    solOpts = optimset('Display','none');
    x1 = fminbnd(fun,lb,ub,solOpts);
    maxQ2 = -fun(x1);
    % 输出
    if silenceOrNot
        return
    end
    msg = sprintf('在TEC片冷热侧温度分别为%.2f[K]和%.2f[K]，输入电流为%.4g[A]时吸热量达最大值%.4g[W]', ...
        Th,Tc,x1,maxQ2);
    logger.trace('TE_ShowPerformance',msg)
    % 绘制不同吸热量时的制冷系数曲线
    current = linspace(0.2*x1,2*x1);
    Q2 = zeros(size(current)); E = zeros(size(current));
    for i = 1:length(current)
        [Q2(i),E(i)] = GetTECHeat(current(i),'cooling',Th,Tc,TEC,TECopts);
    end
    COP2 = Q2./E;
    xAxis = Q2;
    yAxis = COP2;
    plot(xAxis,yAxis)
    xlabel('Q_2 [W]')
    ylabel('COP_2')
    titleStr = sprintf('TEC片冷热侧温度分别为%.2f[K]和%.2f[K]，输入电流：%.4g~%.4g[A]', ...
        Th,Tc,current(1),current(end));
    title(titleStr)