function [QP,STransMembr] = DCMD_SPerm(DirectOpt, MembrProps, SFeedSide, SPermSide)
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
%  QTransMembr - (o real scalar) heat transfer across the membrane
%  STransMembr - (o struct) variable of transmembrane stream
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
lb  = [SPermSide.Temp SPermSide.Temp];
ub  = [SFeedSide.Temp SFeedSide.Temp];
opts = optimoptions(@fmincon,'Display','iter');
[TM,fval,exitflag] = fmincon(fun, TM0, [], [], [], [], lb, ub, [], opts);
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
        STransMembr.Temp = TMH;
    case 1
        STransMembr.Temp = TMC;
end
STransMembr.MassFlow = JM*MembrProps.Area;
STransMembr.MassFraction = 0;
STransMembr.SpecHeat = 4.18e3;
STransMembr.Enthalpy = STransMembr.MassFlow * ...
                       STransMembr.SpecHeat * ...
                       STransMembr.Temp;
QP = JH*MembrProps.Area;
%