%% Define the equation system for coupling DCMD with TEHP
%
%  input arguments
%  T  - (real array) T(1)=TS1H, T(2)=T1H, T(3)=T2C, T(4)=TS2C
%                    T(5)=TS2H, T(6)=T2H, T(7)=T2C, T(8)=TS3C
%                    ...
%  TE - (real array) hot- and cold-side environmental temperatures [K]
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
function [F] = DCMD_EqSys(T, TE, TECs, SInFeeds, SInPerms, Membranes)
%% Check the dimension size of input arguments
% Get the number of DCMD stack
NumStack = length(Membranes);
if (length(SInPerms) == NumStack) && (length(SInFeeds) == NumStack) && ...
   (length(TECs) == NumStack+1) && (length(T)/NumStack == 4)
    InputOK = 1;
else
    InputOK = 0;
end
%% Initialize
fvals = zeros(size(T));
TECold = TE(1);
TEHot = TE(2);
%% Define the equations for energy conservation
if InputOK
    % Energy conservation of TEC1
    fvals(1) = DCMD_EqS(T(1), TECold, TECs(1), 0);
    % Heat balance of feed-side stream 1
    fvals(2) = DCMD_EqF(T(2), T(3), T(1), TECold, TECs(1), SInFeeds(1), ...
                    SInPerms(1), Membranes(1));
    % Heat balance of permeate-side stream 1
    fvals(3) = DCMD_EqP(T(2), T(3), TEHot, T(4), TECs(2), SInFeeds(1), ...
                    SInPerms(1), Membranes(1));
    % Energy conservation of TEC2
    fvals(4) = DCMD_EqS(TEHot, T(4), TECs(2), 0);
end
%
F = fvals;
end