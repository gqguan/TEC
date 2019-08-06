function [ alpha, rho, kappa ] = TE_MaterialProp( T_avg )
%% calculate temperature dependent properties
%   Notes of I/O arguments
%   T_avg - (i double scalar) average temperature [K]
%   alpha - (o double array, length(2)) Seebeck coefficient [V/K]
%   rho   - (o double array, length(2)) electrical conductivity [ohm/m]
%   kappa - (o double array, length(2)) thermal conductivity [W/m-K]
%                            1: p-polar
%                            2: n-polar
%
%   by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
%%  function body
%   parameters of property calculation
params = [2.2224e-5,  9.306e-7, -9.905e-10; ...
           5.112e-7,  1.634e-8,  6.279e-11; ...
           6.2605  , -2.777e-2,  4.131e-5 ];
Ts = [1, T_avg, T_avg^2];
alpha = dot(params(1,:), Ts);
rho   = dot(params(2,:), Ts);
kappa = dot(params(3,:), Ts);
end

