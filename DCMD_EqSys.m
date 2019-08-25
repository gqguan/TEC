%% Define the equation system for coupling DCMD with TEHP
%
%  input arguments
%  TIn- (real array) T(1)=TS1H, T(2)=T1H, T(3)=TM1H, T(4)=TM1C, T(5)=T1C,
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
T = reshape(TIn, [6,NumStack]);
%% Define the equations for energy conservation
for i = 1:NumStack
    % Boundary layer adhered the hot side of TECs(i)
    if i == 1
        Q(i,:) = TE_Heat(T(1,i), TEXs(1), TECs(i));
    else
        Q(i,:) = TE_Heat(T(1,i), T(6,i-1), TECs(i));
    end
    JHSH(i) = Q(i,1)/TECs(i).HTArea;
    [~, hF(i)] = DCMD_TM(SInFeeds(i), -JHSH(i)); % negative JH1H indicates heat flowing into the feed
    fvals(1+(i-1)*6) = DCMD_HeatBalance_BL(T(1,i), T(2,i), JHSH(i), hF(i));
    % Feed-side bulk flow of stage 1
    SOutFeeds(i) = SInFeeds(i); SOutPerms(i) = SInPerms(i);
    SOutFeeds(i).Temp = T(2,i); SOutPerms(i).Temp = T(5,i);
    SOutFeeds(i) = DCMD_PackStream(SOutFeeds(i));
    SOutPerms(i) = DCMD_PackStream(SOutPerms(i));
    [QM(i), SM(i)] = DCMD_Permeation(-1, Membranes(i), SOutFeeds(i), SOutPerms(i));
    fvals(2+(i-1)*6) = DCMD_HeatBalance_BF(T(2,i), Q(i,1), QM(i), SInFeeds(i), SM(i));
    % Boundary layer adhered the hot side of the membrane in the 1st stage
    JHM(i) = QM(i)/Membranes(i).Area;
    fvals(3+(i-1)*6) = DCMD_HeatBalance_BL(T(2,i), T(3,i), JHM(i), hF(i));
    % Boundary layer adhered the cold side of the membrane in the 1st stage
    [QM(i), SM(i)] = DCMD_Permeation(1, Membranes(i), SOutFeeds(i), SOutPerms(i));
    JHM(i) = QM(i)/Membranes(i).Area;
    [~, hP(i)] = DCMD_TM(SInPerms(i), JHM(i));
    fvals(4+(i-1)*6) = DCMD_HeatBalance_BL(T(4,i), T(5,i), JHM(i), hP(i));
    % Permeate-side bulk flow of stage 2
    if i == NumStack
        Q(i+1,:) = TE_Heat(TEXs(2), T(6,i), TECs(i+1));
    else
        Q(i+1,:) = TE_Heat(T(1,i+1), T(6,i), TECs(i+1));
    end
    JHSC(i+1) = Q(i+1,2)/TECs(i+1).HTArea;
    fvals(5+(i-1)*6) = DCMD_HeatBalance_BF(T(5,i), Q(i+1,2), QM(i), SInPerms(i), SM(i));
    % Boundary layer adhered the cold side of TEC2
    fvals(6+(i-1)*6) = DCMD_HeatBalance_BL(T(5,i), T(6,i), JHSC(i+1), hP(i));
end
%     % Boundary layer adhered the hot side of TEC2
%     JHSH(2) = Q(2,1)/TECs(2).HTArea;
%     [~, hF(2)] = DCMD_TM(SInFeeds(2), -JHSH(2)); % negative JH2H indicates heat flowing into the feed
%     fvals(7) = DCMD_HeatBalance_BL(T(7), T(8), JHSH(2), hF(2));
%     % Feed-side bulk flow of stage 2
%     SOutFeeds(2) = SInFeeds(2); SOutPerms(2) = SInPerms(2);
%     SOutFeeds(2).Temp = T(8); SOutPerms(2).Temp = T(11);
%     SOutFeeds(2) = DCMD_PackStream(SOutFeeds(2));
%     SOutPerms(2) = DCMD_PackStream(SOutPerms(2));
%     [QM(2), SM(2)] = DCMD_Permeation(-1, Membranes(2), SOutFeeds(2), SOutPerms(2));
%     fvals(8) = DCMD_HeatBalance_BF(T(8), Q(2,1), QM(2), SInFeeds(2), SM(2));
%     % Boundary layer adhered the hot side of the membrane in the 2nd stage
%     JHM(2) = QM(2)/Membranes(2).Area;
%     fvals(9) = DCMD_HeatBalance_BL(T(8), T(9), JHM(2), hF(2));
%     % Boundary layer adhered the cold side of the membrane in the 2nd stage
%     [QM(2), SM(2)] = DCMD_Permeation(1, Membranes(2), SOutFeeds(2), SOutPerms(2));
%     JHM(2) = QM(2)/Membranes(2).Area;
%     [~, hP(2)] = DCMD_TM(SInPerms(2), JHM(2));
%     fvals(10) = DCMD_HeatBalance_BL(T(10), T(11), JHM(2), hP(2));
%     % Permeate-side bulk flow of stage 2
%     Q(3,:) = TE_Heat(TEXs(2), T(12), TECs(3));
%     JHSC(3) = Q(3,2)/TECs(3).HTArea;
%     fvals(11) = DCMD_HeatBalance_BF(T(11), Q(3,2), QM(2), SInPerms(2), SM(2));
%     % Boundary layer adhered the cold side of TEC3
%     fvals(12) = DCMD_HeatBalance_BL(T(11), T(12), JHSC(3), hP(2));
%% Output results
F = fvals;
end