function [fval] = DCMD_Diff_THTC(x, QH, QC, SInFeed, SInPerm, MembrProps)
%% calculate temperature differences of TH and TC
%  notes of I/O arguments
%  x    - (i real arrays(2)) x(1) = TH, mean temperature of feed-side
%                                       channel in DCMD module [K]
%                            x(2) = TC, mean temperature of permeate-side
%                                       channel in DCMD module [K]
%  QH   - (i real scalar) heat absorbed by feed from the external heat
%                         source [W]
%  QC   - (i real scalar) heat released from the permeate to the external
%                         heat sink [W]
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
%  fval - (o real scalar) differences of [TH TC]
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-15
%  
% initialize
SFeedSide.Temp = x(1);
SPermSide.Temp = x(2);
SFeedSide.MassFraction = SInFeed.MassFraction;
SPermSide.MassFraction = SInPerm.MassFraction;
% feed-side heat and mass balance
DirectOpt = -1;
[QTransMembr,STransMembr] = DCMD_SPerm(DirectOpt, MembrProps, ...
                                       SFeedSide, SPermSide);
SOutFeed = DCMD_SOut(DirectOpt, SInFeed, STransMembr, QH, QTransMembr);
% permeate-side heat and mass balance
DirectOpt = 1;
[QTransMembr,STransMembr] = DCMD_SPerm(DirectOpt, MembrProps, ...
                                       SFeedSide, SPermSide);
SOutPerm = DCMD_SOut(DirectOpt, SInPerm, STransMembr, QTransMembr, QC);
% output
xout = [SOutFeed.Temp,SOutPerm.Temp];
fval = norm(xout-x);
%
end