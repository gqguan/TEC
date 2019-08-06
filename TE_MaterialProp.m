function [ alpha, R, K ] = TE_MaterialProp( T_avg, gf )
%% calculate temperature dependent properties
%   Notes of I/O arguments
%   T_avg - (i double scalar) average temperature [K]
%   gf    - (i double scalar) geometric factor
%   alpha - (o double scalar) Seebeck coefficient [V/K]
%   R     - (o double scalar) electrical conductance [ohm]
%   K     - (o double scalar) thermal conductance [W/K]
%
%   by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
%%  function body
% intialize
results = zeros(0,3);
% parameters of property calculation
params = [2.2224e-5,  9.306e-7, -9.905e-10; ...
           5.112e-7,  1.634e-8,  6.279e-11; ...
           6.2605  , -2.777e-2,  4.131e-5 ];
% temperature array
Ts = [1, T_avg, T_avg^2];
% calculate the properties of thermocouple
for i = 1:3
    results(i) = dot(params(i,:), Ts); % results = [alpha rho kappa]
end
% pack the results
alpha = (results(1)+results(1));
R     = (results(2)+results(2))/gf;
K     = (results(3)+results(3))*gf;
%
end

