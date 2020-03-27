function [T] = TE_JunctionT(Th, Tc, I, TEC)
%% calculate the hot and cold junction temperatures
%  notes of I/O arguments
%  Th  - (i double scalar) hot-side temperature of heat sink [K]
%  Tc  - (i double scalar) cold-side temperature of heat source [K]
%  I   - (i double scalar) input electrical current of TEHP [A]
%  TEC - (i struc) struc variable
%        NumTC     : Number of thermocouples in TEC
%        NumRatio  : ratio of thermocouples in the 1-stage TEC to those in
%                     the 2-stage TEC
%        GeomFactor: geometry factor of thermcouples in TEC [m]
%        HTCoefficient     : Overall heat transfer coefficient [W/m2-K]
%        HTArea            : heat transfer area [m2]
%        SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC
%        ElecConductance   : electrical conductance of 1 and 2 stage of TEC
%        ThermConductance  : thermal conductance of 1 and 2 stage of TEC
%  T   - (o double array(2)) T(1) hot junction temperature [K]
%                            T(2) cold junction temperature [K]
%
%  ## References
%   * Kaushik et al. International Journal of Heat and Mass Transfer
%   86(2015) 843-852
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-11
%
%%
% initialize
syms QH QC T1 T2;
n = TEC.NumTC;
Uh = TEC.HTCoefficient;
Uc = TEC.HTCoefficient;
Ah = TEC.HTArea;
Ac = TEC.HTArea;
[a,R,K] = TE_MaterialProp((Th+Tc)/2, TEC.GeomFactor);
% define the energy balance equations
eq47 = QH == n*(a*I*T1+I^2*R/2-K*(T1-T2));
eq48 = QC == n*(a*I*T2-I^2*R/2-K*(T1-T2));
% substitute eq49 into eq47
eq47a = subs(eq47, QH, Uh*Ah*(T1-Th));
% substitute eq50 into eq48
eq48a = subs(eq48, QC, Uc*Ac*(Tc-T2));
% solve T1 T2
[T1,T2] = solve([eq47a,eq48a], [T1,T2]);
% T1a = ((Uh*Ah*Th+0.5*n*I^2*R)*(Ac*Uc+n*K+n*a*I) + ...
%        (0.5*n^2*K*R*I^2+n*K*Uc*Ac*Tc)) / ...
%       ((n*K+Ac*Uc+n*a*I)*(Ah*Uh+n*K-n*a*I)-n^2*K^2);
% T2a = ((Uc*Ac*Tc+0.5*n*I^2*R)*(Ah*Uh+n*K-n*a*I) + ...
%        (0.5*n^2*K*R*I^2+n*K*Uh*Ah*Th)) / ...
%       ((n*K+Ac*Uc+n*a*I)*(Ah*Uh+n*K-n*a*I)-n^2*K^2);
if T1-T2 > 80
    fprintf('(T1-T2) = %5.2f K is too large!\n', double(T1-T2));
    T = [Th,Tc];
    return
else
    T = eval([T1,T2]);
end
end

