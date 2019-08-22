%% Calculate the difference of TS2H 
%
%  input arguments
%  TS2H      - (real scalar) hot-side temperature of TEC2 [K]
%  TEXs      - (real array) external heat source and sink temperatures [K]
%  TECs      - (struct array) properties of TECs TECs(1) = TEC1
%                                                TECs(2) = TEC2
%                                                TECs(3) = TEC3
%  SInFeeds  - (struct array) properties of feed-side influents
%                       SInFeeds(1): feed-side influent of the 1st stage
%                       SInFeeds(2): feed-side influent of the 2nd stage
%  SInPerms  - (struct array) properties of permeate-side influents
%                       SInPerms(1): permeate-side influent of 1st stage
%                       SInPerms(2): permeate-side influent of 2nd stage
%  Membranes - (struct array) properties of membranes
%                             Membranes(1): membrane of 1st stage
%                             Membranes(2): membrane of 2nd stage
%  output arguments
%  dTS2H     - (real scalar)  
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-21
%
function [dTS2H] = DCMD_Diff_TS2H(TS2H, TEXs, TECs, SInFeeds, SInPerms, ...
                                  Membranes)
%% Initialize
T0 = [SInFeeds(1).Temp+1; SInFeeds(1).Temp; SInPerms(1).Temp; ...
      SInPerms(1).Temp-1];
TE1 = [TEXs(1); TS2H];
%% Calculate the heat flow (QS2H) from 1st stage to 2nd stage
% solve the temperatures in the 1st stage
fun1 = @(T1)DCMD_EqSys(T1, TE1, TECs(1:2), SInFeeds(1), SInPerms(1), ...
                       Membranes(1));
T1 = fsolve(fun1, T0);
% check temperature profile
if (T1(1)>T1(2)) && (T1(2)>T1(3)) && (T1(3)>T1(4))
    fprintf('Temperature profile checked ...\n');
else
    fprintf('Abnormal temperature profile!\n');
end
% solve the temperatures in the 2nd stage
TE2 = [T1(4); TEXs(2)];
fun2 = @(T2)DCMD_EqSys(T2, TE2, TECs(2:3), SInFeeds(2), SInPerms(2), ...
                       Membranes(1));
T2 = fsolve(fun2, T0);
%% Output
dTS2H = T2(1)-TS2H;
%
end