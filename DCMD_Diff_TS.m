function [fval] = DCMD_Diff_TS(x, y, I1, TEC1, I2, TEC2, SInFeed, SInPerm, MembrProps)
%% calculate the temperature differences of TS
%  notes of I/O arguments
%  x - (i real arrays(2)) x(1)=TS1H, T at hot surface of TEC1 [K]
%                         x(2)=TS2C, T at cold surface of TEC2 [K]
%  y - (i real arrays(2)) y(1)=TS1C, T at cold surface of TEC1 [K]
%                         y(2)=TS2H, T at hot surface of TEC2 [K]
%  I1   - (i real scalar) electrical current of TEC1 [A]
%  TEC1 - (i struc) properties of TEC1
%      .NumTC     : Number of thermocouples in TEC1
%      .NumRatio  : ratio of thermocouples in the 1-stage TEC1 to those in
%                   the 2-stage TEC1
%      .GeomFactor: geometry factor of thermcouples in TEC1 [m]
%      .SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC1
%      .ElecConductance   : electrical conductance of 1 and 2 stage of TEC1
%      .ThermConductance  : thermal conductance of 1 and 2 stage of TEC1
%  I2   - (i real scalar) electrical current of TEC2 [A]
%  TEC2 - (i struc) properties of TEC2
%      .NumTC     : Number of thermocouples in TEC2
%      .NumRatio  : ratio of thermocouples in the 1-stage TEC2 to those in
%                   the 2-stage TEC2
%      .GeomFactor: geometry factor of thermcouples in TEC2 [m]
%      .SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC2
%      .ElecConductance   : electrical conductance of 1 and 2 stage of TEC2
%      .ThermConductance  : thermal conductance of 1 and 2 stage of TEC2
%  SInFeed - (i struct) stream properties of feed-side influent
%         .Temp: temperature [K]
%         .MassFlow: mass flowrate [kg/s]
%         .MassFraction: mass fraction of NaCl 
%         .SpecHeat: specific heat [J/kg-K]
%         .Enthalpy: enthalpy [W]
%  SInPerm - (i struct) stream properties of permeate-side influent
%         .Temp: temperature [K]
%         .MassFlow: mass flowrate [kg/s]
%         .MassFraction: mass fraction of NaCl 
%         .SpecHeat: specific heat [J/kg-K]
%         .Enthalpy: enthalpy [W]
%  MembrProps - (i struct) variable of membrane properties
%            .TMH: hot-side temperature of membrane [K]
%            .TMC: cold-side temperature of membrane [K]
%            .Area: effective area of membrane [m2]
%            .Thickness: thickness of membrane [m]
%            .MDCoefficient: MD coefficient [kg/s-m2-Pa]
%            .ThermConductivity: thermal conductivity of membrane
%  fval - (o real scalar) differences of [TS1H TS2C]
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-15
%  
% initialize
TS1H = x(1); TS2C = x(2);
TS1C = y(1); TS2H = y(1);
% calculate the Q1H and Q2C
Q1 = TE_Heat(TS1H, TS1C, I1, TEC1);
Q1H = Q1(1);
Q2 = TE_Heat(TS2H, TS2C, I2, TEC2);
Q2C = Q2(2);
% calculate the T1H and T1C by optimizing feed- and permeate-side 
% temperatures to minimize DCMD_Diff_THTC()
x0 = [SInFeed.Temp SInPerm.Temp];
fun = @(x)DCMD_Diff_THTC(x, Q1H, Q2C, SInFeed, SInPerm, MembrProps);
[T,fval,exitflag] = fminsearch(fun, x0);
T1H = T(1); T1C = T(2);
% calculate the temperatures at both hot surface of TEC1 and cold surface 
% of TEC2
SFeedSide = SInFeed; % UNDER CONSTRUCTION
SPermSide = SInPerm; % UNDER CONSTRUCTION
SFeedSide.Temp = T1H;
SFeedSide = DCMD_PackStream(SFeedSide);
SPermSide.Temp = T1C;
SPermSide = DCMD_PackStream(SPermSide);
xnew(1) = DCMD_TM(SFeedSide, -Q1H/MembrProps.Area);
xnew(2) = DCMD_TM(SPermSide,  Q2C/MembrProps.Area);
fval = norm(x-xnew);
%
end