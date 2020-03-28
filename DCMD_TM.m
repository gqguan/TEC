function [TM, h] = DCMD_TM(Stream, JH)
%% calculate the temperature at the membrane surface
%  notes of I/O arguments
%  input arguments
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
%  output arguments
%  TM     - (o real scalar) wall temperature [K]
%  h      - (o real scalar) overall heat transfer coefficient of boundary
%                           layer [W/m2-K]
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-14
%
%  2019-08-24: add output of heat transfer coefficient
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
% 
% eq = Nu == 3.66+0.19*Gz^0.8/(1+0.117*Gz^0.467); % Guan2012IECR eq.(6)
Re = L*u*rho/mu;
Pr = cp*mu/k;
Gz = Dh/0.006*Re*Pr;
% calculate total heat transfer coefficient
h = k/L*(3.66+0.19*Gz^0.8/(1+0.117*Gz^0.467));
% output
TM = T-JH/h;
%
end