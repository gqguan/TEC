function [fval] = DCMD_Diff_TM(TM, SFeedSide, SPermSide, MembrProps)
%% calculate the differences of [TMH TMC]
%  notes of I/O arguments
%  TM - (i real array(2)) TM(1) = TMH
%                         TM(2) = TMC
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
%  MembrProps - (i struct) variable of membrane properties
%            .TMH: hot-side temperature of membrane [K]
%            .TMC: cold-side temperature of membrane [K]
%            .Area: effective area of membrane [m2]
%            .Thickness: thickness of membrane [m]
%            .MDCoefficient: MD coefficient [kg/s-m2-Pa]
%            .ThermConductivity: thermal conductivity of membrane
%  fval       - (o real scalar) differences of [TMH TMC]
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-15
%
%%
% intialize
TMH = TM(1);
TMC = TM(2);
MF = SFeedSide.MassFraction;
C  = MembrProps.MDCoefficient;
K  = MembrProps.ThermConductivity;
d  = MembrProps.Thickness;
% get permeation flux according to the difference of vapor pressure
PSH = DCMD_SatVapPressure(TMH, MF);
PSC = DCMD_SatVapPressure(TMC);
JM  = C*(PSH-PSC);
% get temperature-dependent latent heat
dHv = DCMD_LatentHeat((TMH+TMC)/2);
% calculate the transmembrane heat flux
JH = JM*dHv+K/d*(TMH-TMC);
TM_New(1) = DCMD_TM(SFeedSide, JH);
TM_New(2) = DCMD_TM(SPermSide, -JH);
fval = norm(TM-TM_New);
end