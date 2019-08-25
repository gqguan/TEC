%% Define the equation system for coupling DCMD with TEHP
%
%  input arguments
%  T  - (real array) T(1)=TS1H, T(2)=T1H, T(3)=TM1H, T(4)=TM1C, T(5)=T1C,
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
%  F  - (real array) left hand side of eq.S, F and P 
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-19
%  
function [F] = DCMD_EqSys(T, TEXs, TECs, SInFeeds, SInPerms, Membranes)
%% Check the dimension size of input arguments
% Get the number of DCMD stack
NumStack = length(Membranes);
if (length(SInPerms) == NumStack) && (length(SInFeeds) == NumStack) && ...
   (length(TECs) == NumStack+1) && (length(T)/NumStack == 6)
    InputOK = 1;
else
    InputOK = 0;
    fprintf('Error of input arguments in DCMD_EqSys()! \n');
end
%% Initialize
fvals = zeros(size(T));
%% Define the equations for energy conservation
if InputOK
    % Boundary layer adhered the hot side of TEC1
    Q1 = TE_Heat(T(1), TEXs(1), TECs(1));
    JH1H = Q1(1)/TECs(1).HTArea;
    [~, hF(1)] = DCMD_TM(SInFeeds(1), -JH1H); % negative JH1H indicates  
                                              % heat flowing into the feed
    fvals(1) = DCMD_HeatBalance_BL(T(1), T(2), JH1H, hF(1));
    % Feed-side bulk flow of stage 1
    SFeedSide = SInFeeds(1); SPermSide = SInPerms(1);
    SFeedSide.Temp = T(2); SPermSide.Temp = T(5);
    SFeedSide = DCMD_PackStream(SFeedSide);
    SPermSide = DCMD_PackStream(SPermSide);
    [QM(1), SM(1)] = DCMD_Permeation(-1, Membranes(1), SFeedSide, SPermSide);
    fvals(2) = DCMD_HeatBalance_BF(T(2), Q1(1), QM(1), SInFeeds(1), SM(1));
    % Boundary layer adhered the hot side of the membrane in the 1st stage
    JH1M = QM(1)/Membranes(1).Area;
    fvals(3) = DCMD_HeatBalance_BL(T(2), T(3), JH1M, hF(1));
    % Boundary layer adhered the cold side of the membrane in the 1st stage
    [QM(1), SM(1)] = DCMD_Permeation(1, Membranes(1), SFeedSide, SPermSide);
    JH1M = QM(1)/Membranes(1).Area;
    [~, hP(1)] = DCMD_TM(SInPerms(1), JH1M);
    fvals(4) = DCMD_HeatBalance_BL(T(4), T(5), JH1M, hP(1));
    % Permeate-side bulk flow of stage 2
    Q2 = TE_Heat(T(7), T(6), TECs(2));
    JH2C = Q2(2)/TECs(2).HTArea;
    fvals(5) = DCMD_HeatBalance_BF(T(5), Q2(2), QM(1), SInPerms(1), SM(1));
    % Boundary layer adhered the cold side of TEC2
    fvals(6) = DCMD_HeatBalance_BL(T(5), T(6), JH2C, hP(1));
    % Boundary layer adhered the hot side of TEC2
    JH2H = Q2(1)/TECs(2).HTArea;
    [~, hF(2)] = DCMD_TM(SInFeeds(2), -JH2H); % negative JH2H indicates  
                                              % heat flowing into the feed
    fvals(7) = DCMD_HeatBalance_BL(T(7), T(8), JH2H, hF(2));
    % Feed-side bulk flow of stage 2
    SFeedSide = SInFeeds(2); SPermSide = SInPerms(2);
    SFeedSide.Temp = T(8); SPermSide.Temp = T(11);
    SFeedSide = DCMD_PackStream(SFeedSide);
    SPermSide = DCMD_PackStream(SPermSide);
    [QM(2), SM(2)] = DCMD_Permeation(-1, Membranes(2), SFeedSide, SPermSide);
    fvals(8) = DCMD_HeatBalance_BF(T(8), Q2(1), QM(2), SInFeeds(2), SM(2));
    % Boundary layer adhered the hot side of the membrane in the 2nd stage
    JH2M = QM(2)/Membranes(2).Area;
    fvals(9) = DCMD_HeatBalance_BL(T(8), T(9), JH2M, hF(2));
    % Boundary layer adhered the cold side of the membrane in the 2nd stage
    [QM(2), SM(2)] = DCMD_Permeation(1, Membranes(2), SFeedSide, SPermSide);
    JH2M = QM(2)/Membranes(2).Area;
    [~, hP(2)] = DCMD_TM(SInPerms(2), JH2M);
    fvals(10) = DCMD_HeatBalance_BL(T(10), T(11), JH2M, hP(2));
    % Permeate-side bulk flow of stage 2
    Q3 = TE_Heat(TEXs(2), T(12), TECs(3));
    JH3C = Q3(2)/TECs(3).HTArea;
    fvals(11) = DCMD_HeatBalance_BF(T(11), Q3(2), QM(2), SInPerms(2), SM(2));
    % Boundary layer adhered the cold side of TEC3
    fvals(12) = DCMD_HeatBalance_BL(T(11), T(12), JH3C, hP(2));
end
%
F = fvals;
end