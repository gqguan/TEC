%% Define the equation system for coupling DCMD with TEHP
%
%  input arguments
%  TIn- (real array(NumStage*6)) 长度为NumStage*6的温度向量，其中每级DCMD包括
%  6个温度（依次分别为热侧壁面、热侧主体、热侧膜面、冷侧膜面、冷侧主体和冷侧壁面温度） 
%                    T(1)=TS1H, T(2)=T1H, T(3)=TM1H, T(4)=TM1C, T(5)=T1C,
%                    T(6)=TS2C, T(7)=TS2H, T(8)=T2H, T(9)=TM2H, 
%                    T(10)=TM2C, T(11)=T2C, T(12) = TS3C
%                    ...
%  TEXs - (real array) hot- and cold-side environmental temperatures [K]
%  TECs - (struct array) properties of TECs TECs(1) = TEC1
%                                           TECs(2) = TEC2
%                                           TECs(3) = TEC3
%                                           ...
%  SInFeeds - (struct array) properties of feed-side influents
%                       SInFeeds(i): feed-side influent of i-th stage
%  SInPerms - (struct array) properties of permeate-side influents
%                       SInFeeds(i): permeate-side influent of i-th stage
%  Membranes - (struct array) properties of membranes
%                             Membranes(i): membrane of i-th stage
%  output arguments
%  F  - (real array(12)) left hand side of eqs.(1)-(12) 
%  Q  - (real matrix(3,2)) absorbed and released heats of TEC1-3
%  QM - (real array(2)) transmembrane heats
%  SM - (struct array(2)) properties of permeation streams
%  SOutFeeds - (struct array(2)) properties of feed-side effluents
%  SOutPerms - (struct array(2)) properties of permeate-side effluents
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-23
%  
function [F, Q, QM, SM, SOutFeeds, SOutPerms] = ...
         DCMD_EqSys(TIn, TEXs, TECs, SInFeeds, SInPerms, Membranes)
%% Check the dimension size of input arguments
% Get the number of DCMD stack
NumStack = length(Membranes);
if (length(SInPerms) == NumStack) && (length(SInFeeds) == NumStack) && ...
   (length(TECs) == NumStack+1) && (length(TIn)/NumStack == 6)
    InputOK = 1;
else
    InputOK = 0;
    fprintf('Error of input arguments in DCMD_EqSys()! \n');
    return;
end
%% Initialize
fvals     = zeros(size(TIn));
Q         = zeros(length(TECs),2);
JHSH      = zeros(size(TECs));
JHSC      = zeros(size(TECs));
JHM       = zeros(size(Membranes));
hF        = zeros(size(Membranes));
hP        = zeros(size(Membranes));
QM        = zeros(size(Membranes));
SOutFeeds = struct(SInFeeds);
SOutPerms = struct(SInPerms);
SM        = struct(SInPerms);
T = reshape(TIn, [6,NumStack]); % 将输入的温度向量重整为尺寸是[6,NumStack]的二维数组
%% Define the equations for energy conservation
for i = 1:NumStack
    % Boundary layer adhered the hot side of TECs(i)
    % 计算该级膜组件冷热侧TEC的吸放热量
    if i == 1
        Q(i,:) = TE_Heat(T(1,i), TEXs(1), TECs(i));
    else
        Q(i,:) = TE_Heat(T(1,i), T(6,i-1), TECs(i));
    end
    JHSH(i) = Q(i,1)/TECs(i).HTArea; % 从壁面传入热侧流体的热通量
    % 计算热侧壁面与主流间的传热系数（W/m2-K）
    [~, hF(i)] = DCMD_TM(SInFeeds(i), -JHSH(i)); % negative JH1H indicates heat flowing into the feed
    % 计算热侧壁面传热通量的差值，用于求解热侧壁面温度使壁面向主流传热通量与壁面流出热量相等
    fvals(1+(i-1)*6) = DCMD_HeatBalance_BL(T(1,i), T(2,i), JHSH(i), hF(i))*Membranes(i).Area; % 注意该函数计算结果为热通量差值
    % Feed-side bulk flow of stage i
    % 初始化该级DCMD膜组件流出的冷热侧物流
    SOutFeeds(i) = SInFeeds(i); SOutPerms(i) = SInPerms(i);
    % 设定冷热侧输出物流的温度
    SOutFeeds(i).Temp = T(2,i); SOutPerms(i).Temp = T(5,i);
    % 重新计算物流性质
    SOutFeeds(i) = DCMD_PackStream(SOutFeeds(i));
    SOutPerms(i) = DCMD_PackStream(SOutPerms(i));
    % 设定膜两侧温度
    Membranes(i).TMH = T(3,i);
    Membranes(i).TMC = T(4,i);
    % 计算热侧流向膜面的热流和物流量
    [QM(i), SM(i)] = DCMD_Permeation(-1, Membranes(i), SOutFeeds(i), SOutPerms(i), 1);
    % 计算热侧主流CV中的能量差值，用于求解热侧主体温度使进出物流的焓变等于CV的净输入热量
    fvals(2+(i-1)*6) = DCMD_HeatBalance_BF(T(2,i), Q(i,1), QM(i), SInFeeds(i), SM(i));
    % Boundary layer adhered the hot side of the membrane in the i-th stage
    % 计算热侧从主流向膜面的传热通量
    JHM(i) = QM(i)/Membranes(i).Area;
    % 计算热侧主流与膜面传热通量的差值，用于求解热侧膜面温度使主流向膜面边界层流出的热量等于膜面边界层内的传热量
    fvals(3+(i-1)*6) = DCMD_HeatBalance_BL(T(2,i), T(3,i), JHM(i), hF(i))*Membranes(i).Area;
    % Boundary layer adhered the cold side of the membrane in the 1st stage
    % 计算冷侧从膜面流向主流的热流和物流量
    [QM(i), SM(i)] = DCMD_Permeation(1, Membranes(i), SOutFeeds(i), SOutPerms(i), 1);
    % 计算冷侧从膜面向主流的传热通量
    JHM(i) = QM(i)/Membranes(i).Area;
    % 计算冷侧膜面边界层的传热系数
    [~, hP(i)] = DCMD_TM(SInPerms(i), JHM(i));
    % 计算冷侧边界层中的净热量差，用于求解冷侧膜面温度使膜面向主流的传热量等于跨膜传热量
    fvals(4+(i-1)*6) = DCMD_HeatBalance_BL(T(4,i), T(5,i), JHM(i), hP(i))*Membranes(i).Area;
    % Permeate-side bulk flow of stage i+1
    % 计算第i级DCMD中冷侧TEC（即第i+1级DCMD中的热侧TEC）壁面的吸放热量
    if i == NumStack
        Q(i+1,:) = TE_Heat(TEXs(2), T(6,i), TECs(i+1));
    else
        Q(i+1,:) = TE_Heat(T(1,i+1), T(6,i), TECs(i+1));
    end
    % 计算第i级DCMD中冷侧壁面的吸热量
    JHSC(i+1) = Q(i+1,2)/TECs(i+1).HTArea;
    % 计算冷侧主流CV的净热量差，用于求解冷侧主流温度使进出物流的焓变等于CV的净输入热量
    fvals(5+(i-1)*6) = DCMD_HeatBalance_BF(T(5,i), Q(i+1,2), QM(i), SInPerms(i), SM(i));
    % Boundary layer adhered the cold side of TEC2
    % 计算冷侧壁面边界层的净热量差，用于求解冷侧壁面温度使边界层内传热量等于第i+1级DCMD热侧TEC的冷侧吸热量
    fvals(6+(i-1)*6) = DCMD_HeatBalance_BL(T(5,i), T(6,i), JHSC(i+1), hP(i))*Membranes(i).Area;
end

%% Output results
F = fvals;
end