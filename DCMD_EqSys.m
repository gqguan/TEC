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
         DCMD_EqSys(T, TEXs, TECs, SInFeeds, SInPerms, Membranes)
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
    Q(1,:) = TE_Heat(T(1), TEXs(1), TECs(1));
    JHSH(1) = Q(1,1)/TECs(1).HTArea;
    [~, hF(1)] = DCMD_TM(SInFeeds(1), -JHSH(1)); % negative JH1H indicates heat flowing into the feed
    fvals(1) = DCMD_HeatBalance_BL(T(1), T(2), JHSH(1), hF(1));
    % Feed-side bulk flow of stage 1
    SOutFeeds(1) = SInFeeds(1); SOutPerms(1) = SInPerms(1);
    SOutFeeds(1).Temp = T(2); SOutPerms(1).Temp = T(5);
    SOutFeeds(1) = DCMD_PackStream(SOutFeeds(1));
    SOutPerms(1) = DCMD_PackStream(SOutPerms(1));
    [QM(1), SM(1)] = DCMD_Permeation(-1, Membranes(1), SOutFeeds(1), SOutPerms(1));
    fvals(2) = DCMD_HeatBalance_BF(T(2), Q(1,1), QM(1), SInFeeds(1), SM(1));
    % Boundary layer adhered the hot side of the membrane in the 1st stage
    JHM(1) = QM(1)/Membranes(1).Area;
    fvals(3) = DCMD_HeatBalance_BL(T(2), T(3), JHM(1), hF(1));
    % Boundary layer adhered the cold side of the membrane in the 1st stage
    [QM(1), SM(1)] = DCMD_Permeation(1, Membranes(1), SOutFeeds(1), SOutPerms(1));
    JHM(1) = QM(1)/Membranes(1).Area;
    [~, hP(1)] = DCMD_TM(SInPerms(1), JHM(1));
    fvals(4) = DCMD_HeatBalance_BL(T(4), T(5), JHM(1), hP(1));
    % Permeate-side bulk flow of stage 2
    Q(2,:) = TE_Heat(T(7), T(6), TECs(2));
    JHSC(2) = Q(2,2)/TECs(2).HTArea;
    fvals(5) = DCMD_HeatBalance_BF(T(5), Q(2,2), QM(1), SInPerms(1), SM(1));
    % Boundary layer adhered the cold side of TEC2
    fvals(6) = DCMD_HeatBalance_BL(T(5), T(6), JHSC(2), hP(1));
    % Boundary layer adhered the hot side of TEC2
    JHSH(2) = Q(2,1)/TECs(2).HTArea;
    [~, hF(2)] = DCMD_TM(SInFeeds(2), -JHSH(2)); % negative JH2H indicates heat flowing into the feed
    fvals(7) = DCMD_HeatBalance_BL(T(7), T(8), JHSH(2), hF(2));
    % Feed-side bulk flow of stage 2
    SOutFeeds(2) = SInFeeds(2); SOutPerms(2) = SInPerms(2);
    SOutFeeds(2).Temp = T(8); SOutPerms(2).Temp = T(11);
    SOutFeeds(2) = DCMD_PackStream(SOutFeeds(2));
    SOutPerms(2) = DCMD_PackStream(SOutPerms(2));
    [QM(2), SM(2)] = DCMD_Permeation(-1, Membranes(2), SOutFeeds(2), SOutPerms(2));
    fvals(8) = DCMD_HeatBalance_BF(T(8), Q(2,1), QM(2), SInFeeds(2), SM(2));
    % Boundary layer adhered the hot side of the membrane in the 2nd stage
    JHM(2) = QM(2)/Membranes(2).Area;
    fvals(9) = DCMD_HeatBalance_BL(T(8), T(9), JHM(2), hF(2));
    % Boundary layer adhered the cold side of the membrane in the 2nd stage
    [QM(2), SM(2)] = DCMD_Permeation(1, Membranes(2), SOutFeeds(2), SOutPerms(2));
    JHM(2) = QM(2)/Membranes(2).Area;
    [~, hP(2)] = DCMD_TM(SInPerms(2), JHM(2));
    fvals(10) = DCMD_HeatBalance_BL(T(10), T(11), JHM(2), hP(2));
    % Permeate-side bulk flow of stage 2
    Q(3,:) = TE_Heat(TEXs(2), T(12), TECs(3));
    JHSC(3) = Q(3,2)/TECs(3).HTArea;
    fvals(11) = DCMD_HeatBalance_BF(T(11), Q(3,2), QM(2), SInPerms(2), SM(2));
    % Boundary layer adhered the cold side of TEC3
    fvals(12) = DCMD_HeatBalance_BL(T(11), T(12), JHSC(3), hP(2));
end
%% Output results
F = fvals;
end