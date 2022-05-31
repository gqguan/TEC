%% 一维集成热泵DCMD膜组件温度分布
%
% by Dr. Guan Guoqiang @ SCUT on 2022-05-02

function [profile,sOut] = TEHPiDCMD(sIn,TECs,TEXs,membrane,flowPattern,opts)
%% 初始化
% 调试用输入参数
if nargin == 0
    [~,s1,membrane] = InitStruct();
    T1 = 273.15+50; T2 = 273.15+45; % [K]
    W1 = 1.217e-2; W2 = 6.389e-3; % [kg/s]
    s1.Temp = T1;
    s1.MassFlow = W1;
    sIn(1) = DCMD_PackStream(s1); % 膜组件料液侧进料
    s2 = s1;
    s2.Temp = T2;
    s2.MassFlow = W2;
    sIn(2) = DCMD_PackStream(s2); % 膜组件渗透侧进料
    iTECs = [20,1];
    T0 = 298.15;
    TEXs = [T0,T0];
    % set properties for all TECs
    load('TEC_Params.mat','TEC_Params') % 载入已有的TEC计算参数
    TECs(1) = TEC_Params.TEC(iTECs(1));
    % TEC计算参数
    opts = [TEC_Params.opt1(iTECs(1)),TEC_Params.opt2(iTECs(1))];
    TECs(2) = TEC_Params.TEC(iTECs(2)); 
    flowPattern = "countercurrent";
    opts = [0 0];
end
nCV = 20;
eps = 1e-8;
% 添加\Common目录以便调用自定义公用函数
homePath = cd;
idxPath = strfind(homePath,'\');
addpath([homePath(1:idxPath(2)),'Common\'])
% 日志文件
logFilename = [homePath,'\',datestr(datetime,'yyyymmdd'),'.log'];
logger = log4m.getLogger(logFilename);
logger.setFilename(logFilename);

%% 检查输入参数
% 检查输入变量数量
if nargin == 6 || nargin == 0
    % 检查输入变量类型
    [flag,msg] = ChkArgType({'struct','struct','double','struct','string','double'}, ...
        sIn,TECs,TEXs,membrane,flowPattern,opts);
    if flag
%         logger.trace('TEHPiDCMD',msg)
    else
        logger.error('TEHPiDCMD',msg)
        return
    end
    % 检查输入变量尺寸
    lenArgs = [2,2,2,1,1,2];
    idxArg = ([length(sIn),length(TECs),length(TEXs),length(membrane), ...
        length(flowPattern),length(opts)] ~= lenArgs);
    if any(idxArg)
        prompt = sprintf('第%d个输入参数尺寸应为%d',find(idxArg,1),lenArgs(find(idxArg,1)));
        logger.error('TEHPiDCMD',prompt)
        return
    end
else
    logger.error('TEHPiDCMD','输入参数数量应为5')
    return
end

%% 

%% Solve temperatures
% 初值
T1 = sIn(1).Temp;
T2 = sIn(2).Temp;
T0 = [T1+1; T1; T1-1; T2+1; T2; T2-1];
% 温度约束条件
lb = ones(size(T0))*273.16;
ub = ones(size(T0))*(273.15+98);
% 求解参数
solOpts = optimoptions(@lsqnonlin, 'Display', 'none');
% 沿流向分M个控制体CV
membr(1:nCV) = membrane;
for iCV = 1:nCV
    membr(iCV).Area = membrane.Area/nCV;
end
localS1(1:nCV+1) = sIn(1); localS2(1:nCV+1) = sIn(2);
localQTEC = cell(1,nCV);
localQM = zeros(1,nCV);
localSM(1:nCV) = sIn(1);
localT = zeros(length(T0),nCV);

switch flowPattern
    case("cocurrent")
        logger.trace('TEHPiDCMD','并流计算')
        for iCV = 1:nCV
            fun = @(T)DCMD_EqSys(T, TEXs, TECs, localS1(iCV), localS2(iCV), membr(iCV), opts);
            T = lsqnonlin(fun, T0, lb, ub, solOpts);
            [~, localQTEC{iCV}, localQM(iCV), localSM(iCV), localS1(iCV+1), localS2(iCV+1)] = fun(T);
            localT(:,iCV) = T;
            if any(ub-T < eps)
                msg = sprintf('料液侧第%d个CV温度过高',iCV);
                logger.warn('TEHPiDCMD',msg)
            end
            if any(T-lb < eps)
                msg = sprintf('料液侧第%d个CV温度过低',iCV);
                logger.warn('TEHPiDCMD',msg)
            end
        end
    case("countercurrent")
        localS2 = fliplr(localS2);
        localS2New = localS2;
        interativeOPT = true;
        nInteration = 1;
        while interativeOPT
            for iCV = 1:nCV
                fun = @(T)DCMD_EqSys(T, TEXs, TECs, localS1(iCV), localS2(iCV+1), membr(iCV), opts);
                T = lsqnonlin(fun, T0, lb, ub, solOpts);
                [~, localQTEC{iCV}, localQM(iCV), localSM(iCV), localS1(iCV+1), localS2New(iCV)] = fun(T);
                localT(:,iCV) = T;
            end
            relDeltaS2T = abs([localS2.Temp]-[localS2New.Temp])./[localS2New.Temp];
            maxRelDT = max(relDeltaS2T);
            if maxRelDT > eps 
                if ChkDivergency(maxRelDT,5)
                    interativeOPT = false;
                    msg = sprintf('迭代求解逆流过程时温度发散（dRelT=%.4e）',maxRelDT);
                    logger.error('TEHPiDCMD',msg)
                end
                localS2 = localS2New;
                nInteration = nInteration+1;
            else
                interativeOPT = false;
%                 msg = sprintf('逆流计算迭代次数%d',nInteration);
%                 logger.trace('TEHPiDCMD',msg)
            end
        end
    otherwise
        logger.error('TEHPiDCMD','第5个输出参数（%s）无效',flowPattern)
end

sOut(1) = localS1(end);
sOut(2) = localS2(end);
profile.T = localT;
profile.QTEC = localQTEC;
profile.QM = localQM;
profile.SM = localSM;
profile.S1 = localS1;
profile.S2 = localS2;
profile.Remarks = sprintf('DCMD膜组件流动模式：%s',flowPattern);

