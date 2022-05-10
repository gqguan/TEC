function out = TE_ShowPerformance(Th,Tc,TEC,TECopts,opStr)
    % 日志文件
    logFilename = [cd,'\',datestr(datetime,'yyyymmdd'),'.log'];
    logger = log4m.getLogger(logFilename);
    logger.setFilename(logFilename);
    % 初始化
    if ~exist('opStr','var')
        opStr = 'silence';
    end
    showFig = false;
    % 函数调试用
    if nargin == 0
        Th = 273.15+45; % [K]
        Tc = 273.15+5; % [K]
        load('TEC_Params.mat', 'TEC_Params')
        TEC = TEC_Params.TEC(14,:);
        TECopts = [1,0];
        opStr = 'max.COP2';
        showFig = true;
    end
    strIU = {'Current','Voltage'};
    strUnit = {'A','V'};
    % 设定电流范围
    lb = 0.5; ub = 10; % 电流范围0.5~10[A]
    % 计算最大吸热量
    fun1 = @(x1)-GetTECHeat(x1,'cooling',Th,Tc,TEC,TECopts); % 负值求最小即求最大吸热量
    solOpts = optimset('Display','none');
    solX1 = fminbnd(fun1,lb,ub,solOpts);
    
    % 计算最大制冷系数
    fun2 = @(x2)-MaxCOP(x2,'cooling',Th,Tc,TEC,TECopts);
    solX2 = fminbnd(fun2,0.12*solX1,2*solX1,solOpts);
    % 输出
    msg = sprintf('在TEC片冷热侧温度分别为%.2f[K]和%.2f[K]',Th,Tc);
    msg = sprintf('%s，输入%s为',msg,strIU{TECopts(2)+1});
    switch opStr 
        case 'silence'
            out = [];
            return
        case 'max.Q2'
            out = -fun1(solX1);
            msg = sprintf('%s%.4g[%s]时吸热量达最大值%.4g[W]', ...
                msg,solX1,strUnit{TECopts(2)+1},out);
            logger.trace('TE_ShowPerformance',msg)
        case 'max.COP2'
            out = -fun2(solX2);
            msg = sprintf('%s%.4g[%s]时制冷系数达最大值%.4g', ...
                msg,solX2,strUnit{TECopts(2)+1},out);
            logger.trace('TE_ShowPerformance',msg)
        case 'max.All'
            out.MaxQ2.Value = -fun1(solX1);
            out.MaxQ2.(strIU{TECopts(2)+1}) = solX1;
            out.MaxCOP2.Value = -fun2(solX2);
            out.MaxCOP2.(strIU{TECopts(2)+1}) = solX2;
            out.MaxCOP2.Q2 = -fun1(solX2);
        otherwise
            error('TE_ShowPerformance()输入操作参数有误！')
    end
    
    if showFig 
        % 绘制不同吸热量时的制冷系数曲线
        current = linspace(0.12*solX1,2*solX1);
        Q2 = zeros(size(current)); E = zeros(size(current));
        for i = 1:length(current)
            [Q2(i),E(i)] = GetTECHeat(current(i),'cooling',Th,Tc,TEC,TECopts);
        end
        COP2 = Q2./E;
        xAxis = Q2;
        yAxis = COP2;
        figure(1)
        plot(xAxis,yAxis)
        xlabel('Q_2 [W]')
        ylabel('COP_2')
        titleStr = sprintf('TEC片冷热侧温度分别为%.2f[K]和%.2f[K]，输入电流：%.4g~%.4g[A]', ...
            Th,Tc,current(1),current(end));
        title(titleStr)
    end
    
    function COP = MaxCOP(inArg,opStr,Th,Tc,TEC,TECopts)
        [heat,power] = GetTECHeat(inArg,opStr,Th,Tc,TEC,TECopts);
        COP = heat/power;
    end
end