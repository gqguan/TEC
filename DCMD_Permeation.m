function [QM,SM] = DCMD_Permeation(DirectOpt, MembrProps, SFeedSide, SPermSide)
%% calculate transmembrane permeation
%  notes of I/O arguments
%  DirectOpt  - (i integer scalar) opt = -1: flow into membrane
%                                         1: flow out from membrane
%  MembrProps - (i struct) variable of membrane properties
%            .TMH: hot-side temperature of membrane [K]
%            .TMC: cold-side temperature of membrane [K]
%            .Area: effective area of membrane [m2]
%            .Thickness: thickness of membrane [m]
%            .MDCoefficient: MD coefficient [kg/s-m2-Pa]
%            .ThermConductivity: thermal conductivity of membrane
%  SFeedSide  - (i struct) properties of feed-side stream
%           .Temp: temperature [K]
%           .MassFlow: mass flowrate [kg/s]
%           .Velocity: velocity [m/s]
%           .MassFraction: mass fraction of NaCl
%           .Density: density [kg/m3]
%           .Viscosity: dynamic viscosity [Pa-s]
%           .SpecHeat: specific heat [J/kg-K]
%           .ThermCond: thermal conductivity [W/m-K]
%           .Enthalpy: enthalpy [W]
%  SPermSide  - (i struct) properties of permeate-side stream
%           .Temp: temperature [K]
%           .MassFlow: mass flowrate [kg/s]
%           .Velocity: velocity [m/s]
%           .MassFraction: mass fraction of NaCl
%           .Density: density [kg/m3]
%           .Viscosity: dynamic viscosity [Pa-s]
%           .SpecHeat: specific heat [J/kg-K]
%           .ThermCond: thermal conductivity [W/m-K]
%           .Enthalpy: enthalpy [W]
%  QM - (o real scalar) heat transfer across the membrane
%  SM - (o struct) variable of transmembrane stream
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-13
%
%%
% intialize
MF = SFeedSide.MassFraction;
C  = MembrProps.MDCoefficient;
K  = MembrProps.ThermConductivity;
d  = MembrProps.Thickness;
% calculate the temperature at the membrane surface 
fun = @(TM)DCMD_Diff_TM(TM, SFeedSide, SPermSide, MembrProps);
TM0 = [SFeedSide.Temp-1 SPermSide.Temp+1];
% % use fmincon() to get the TMs
% lb  = [SPermSide.Temp SPermSide.Temp];
% ub  = [SFeedSide.Temp SFeedSide.Temp];
% opts = optimoptions(@fmincon,'Display','off');
% [TM,fval,exitflag] = fmincon(fun, TM0, [], [], [], [], lb, ub, [], opts);
% use fminsearch() to get the TMs
TM = fminsearch(fun, TM0);
TMH = TM(1); TMC = TM(2);
% get permeation flux according to the difference of vapor pressure
PSH = DCMD_SatVapPressure(TMH, MF);
PSC = DCMD_SatVapPressure(TMC);
JM  = C*(PSH-PSC);
% get temperature-dependent latent heat
dHv = DCMD_LatentHeat((TMH+TMC)/2);
% calculate the transmembrane heat flux
JH = JM*dHv+K/d*(TMH-TMC);
% output
switch DirectOpt
    case -1
        SM.Temp = TMH;
    case 1
        SM.Temp = TMC;
end
SM.MassFlow = DirectOpt*JM*MembrProps.Area;
SM.MassFraction = 0;
SM.SpecHeat = 4.18e3;
SM = DCMD_PackStream(SM);
QM = JH*MembrProps.Area;
%