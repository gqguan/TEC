function [Tm, TEC] = TE_Tm(Th, Tc, I, TEC)
%% calculate the junction temperature between the 2 stages of TEC
%  notes of I/O arguments
%  Th  - (i double scalar) hot-side temperature [K]
%  Tc  - (i double scalar) cold-side temperature [K]
%  I   - (i double scalar) serial electric current flowing through the two
%                         stages [A]
%  TEC - (i struc) struc variable
%       NumTC     : Number of thermocouples in TEC
%       NumRatio  : ratio of thermocouples in the 1-stage TEC to those in
%                     the 2-stage TEC
%       GeomFactor: geometry factor of thermcouples in TEC [m]
%       SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC
%       ElecConductance   : electrical conductance of 1 and 2 stage of TEC
%       ThermConductance  : thermal conductance of 1 and 2 stage of TEC
%  Tm  - (o double scalar) junction temperature [K]
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-07
%  
%  2019-08-10: use fminbnd() to calculate the Tm
%
% initialize
syms Tm;
% use temperature-independant properties at T = (Th+Tc)/2
[a, R, K] = TE_MaterialProp((Th+Tc)/2, TEC.GeomFactor);
a1 = a; a2 = a; R1 = R; R2 = R; K1 = K; K2 = K;
r = TEC.NumRatio;
% calculate the thermocouple number in the first stage of 2-stage TEC
N0 = TEC.NumTC/(TEC.NumRatio+1);
% solve eq.(3)
eq3 = (I*a1*Tm-I^2*R1/2-K1*(Th-Tm))*r == I*a2*Tm+I^2*R2/2-K2*(Tm-Tc);
Tm = double(solve(eq3, Tm));
% % calculate a1 R1 K1
% [a1, R1, K1] = TE_MaterialProp((Th+Tm)/2, TEC.GeomFactor);
% % calculate a2 R2 K2
% [a2, R2, K2] = TE_MaterialProp((Tc+Tm)/2, TEC.GeomFactor);
% 输出热电偶参数
TEC.SeebeckCoefficient = [a1, a2];
TEC.ElecConductance    = [R1, R2];
TEC.ThermConductance   = [K1, K2];
%
end
