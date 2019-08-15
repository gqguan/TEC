function [TM] = DCMD_TM(Stream, JH)
%% calculate the temperature at the membrane surface
%  notes of I/O arguments
%  stream - (i struct) properties of stream
%           .Temp: temperature [K]
%           .MassFlow: mass flowrate [kg/s]
%           .Velocity: velocity [m/s]
%           .MassFraction: mass fraction of NaCl
%           .Density: density [kg/m3]
%           .Viscosity: dynamic viscosity [Pa-s]
%           .SpecHeat: specific heat [J/kg-K]
%           .Enthalpy: enthalpy [W]
%  JH     - (i real scalar) heat flux [W/m2]
%  TM     - (o real scalar)
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-14
%  
% initialize
L   = 0.008; % characteristic length [m]
T   = Stream.Temp;
k   = Stream.ThermCond;
u   = Stream.Velocity;
rho = Stream.Density;
mu  = Stream.Viscosity;
cp  = Stream.SpecHeat;
Dh  = 4*40*6/(2*(40+6))*1e-3; % hydraulic diameter [m]
syms h; % overall heat transfer coefficient [W/m2-K]
syms Nu Re Pr Gz; % dimensionless numbers
% solve h
% eq = Nu == 0.16*Re^0.67*Pr^0.33; % Guan2012IECR eq.(6)
% eq = Nu == 1.3*Re^0.645*Pr^0.38; % Ibrahim2013AIChE tab.1 ref.[29]
eq = Nu == 3.66+0.19*Gz^0.8/(1+0.117*Gz^0.467); % Guan2012IECR eq.(6)
% eq = Nu == 0.74*Re^0.2*(Gz*Pr)^0.1*Pr^0.2; % Tomaszewska2000JMS eq.(13)
eq = subs(eq, Gz, Dh/0.006*Re*Pr);
eq = subs(eq, [Nu,Re,Pr], [h*L/k,L*u*rho/mu,cp*mu/k]);
h = double(solve(eq, h));
% output
TM = T-JH/h;
%
end