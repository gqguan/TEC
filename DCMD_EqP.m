%% Define the energy balance equation for permeate-side stream
%
%  notes of I/O arguments
%  TH  - (i real scalar) bulk temperature of feed-side stream [K]
%  TC  - (i real scalar) bulk temperature of permeate-side stream [K]
%  TSH - (i real scalar) hot-side temperature of TEC [K]
%  TSC - (i real scalar) cold-side temperature of TEC [K]
%  I   - (i real scalar) input electrical current of TEC [A]
%  TEC - (i struc) properties of TEC
%     .NumTC     : Number of thermocouples in TEC
%     .NumRatio  : ratio of thermocouples in the 1-stage TEC to those in
%                  the 2-stage TEC
%     .GeomFactor: geometry factor of thermcouples in TEC [m]
%     .SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC
%     .ElecConductance   : electrical conductance of 1 and 2 stage of TEC
%     .ThermConductance  : thermal conductance of 1 and 2 stage of TEC
%     .Voltage           : electrical voltage [V]
%     .Current           : electrical current [A]
%  SInFeed - (i struc) properties of feed-side Influent
%         .Temp: temperature [K]
%         .MassFlow: mass flowrate [kg/s]
%         .Velocity: velocity [m/s]
%         .MassFraction: mass fraction of NaCl
%         .Density: density [kg/m3]
%         .Viscosity: dynamic viscosity [Pa-s]
%         .SpecHeat: specific heat [J/kg-K]
%         .ThermCond: thermal conductivity [W/m-K]
%         .Enthalpy: enthalpy [W]
%  SInPerm - (i struc) properties of permeate-side Influent
%         .Temp: temperature [K]
%         .MassFlow: mass flowrate [kg/s]
%         .Velocity: velocity [m/s]
%         .MassFraction: mass fraction of NaCl
%         .Density: density [kg/m3]
%         .Viscosity: dynamic viscosity [Pa-s]
%         .SpecHeat: specific heat [J/kg-K]
%         .ThermCond: thermal conductivity [W/m-K]
%         .Enthalpy: enthalpy [W]
%  MembrProps - (i struct) properties of membrane
%            .TMH: hot-side temperature of membrane [K]
%            .TMC: cold-side temperature of membrane [K]
%            .Area: effective area of membrane [m2]
%            .Thickness: thickness of membrane [m]
%            .MDCoefficient: MD coefficient [kg/s-m2-Pa]
%            .ThermConductivity: thermal conductivity of membrane
%  dE  - (o real scalar) differences of energy balance
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-19
%  
function [dE, SOutPerm] = DCMD_EqP(TH, TC, TSH, TSC, TEC, SInFeed, SInPerm, MembrProps)
%% Calculate the heats flowing into the peremate
Q = TE_Heat(TSH, TSC, TEC);
QC = Q(2);
%% Calculate transmembrane flux
DirectOpt = 1;
% set the influent as the bulk-flow stream
SFeedSide = SInFeed;
SFeedSide.Temp = TH;
SFeedSide = DCMD_PackStream(SFeedSide);
SPermSide = SInPerm;
SPermSide.Temp = TC;
SPermSide = DCMD_PackStream(SPermSide);
[QM, SM] = DCMD_SPerm(DirectOpt, MembrProps, SFeedSide, SPermSide);
%% calculate the effluent
SOutPerm = SInPerm;
SOutPerm.Temp = TC;
SOutPerm.MassFlow = SInPerm.MassFlow + SM.MassFlow;
SOutPerm = DCMD_PackStream(SOutPerm);
HIn = SInPerm.Enthalpy;
HOut = SOutPerm.Enthalpy;
HMIn = SM.Enthalpy;
%% conduct the energy balance
dE = HIn+QM+HMIn-QC-HOut;
%
end