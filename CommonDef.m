%% 程序公用变量定义
%
% by Dr. Guan Guoqiang @ SCUT on 2020/05/02
%

%%
clear
% 定义模组件流道几何尺寸
%  DuctGeom - geometric parameters of flowing channel in MD module
%           .Length (real array(2)) length along the flowing direction in 
%                                   both sides of MD module [m]
%           .Height (real array(2)) height of rectanglarly wetted perimeter
%                                   in both sides of MD module [m]
%           .Width (real array(2))  width of rectanglarly wetted perimeter
%                                   in both sides of MD module [m]
DuctGeom = struct('Length', 0.04,  ...
                 'Height', 0.006, ...
                 'Width',  0.04 );
% 定义物料数据结构
Stream = struct('Temp', 295.15,  ...
                'MassFlow', 0.015,  ...
                'MassFraction', 0,  ...
                'Density', 1e3,     ...
                'Viscosity', 1e-3,  ...
                'SpecHeat', 4.18e3, ...
                'ThermCond', 0.6);
% 定义膜材料性质
%  MembrProps.TMH: hot-side temperature of membrane [K]
%            .TMC: cold-side temperature of membrane [K]
%            .Area: effective area of membrane [m2]
%            .Thickness: thickness of membrane [m]
%            .MDCoefficient: MD coefficient [kg/s-m2-Pa]
%            .ThermConductivity: thermal conductivity of membrane
MembrProps = struct('TMH', [], 'TMC', [], 'Area', 0.0016, ...
                    'Thickness', 1.5e-4, 'MDCoefficient', 3.2e-7, ...
                    'ThermConductivity', (0.18*0.3+0.025*0.7));