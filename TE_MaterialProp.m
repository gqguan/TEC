function [ TEC ] = TE_MaterialProp( T_avg, TEC, opt )
%% calculate temperature dependent properties
%  Notes of I/O arguments
%  T_avg - (i double scalar) average temperature [K]
%  TEC   - (i/o struc) struc variable
%     .NumTC     : Number of thermocouples in TEC
%     .NumRatio  : ratio of thermocouples in the 1-stage TEC to those in
%                    the 2-stage TEC
%     .GeomFactor: geometry factor of thermcouples in TEC [m]
%     .SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC
%     .ElecConductance   : electrical conductance of 1 and 2 stage of TEC
%     .ThermConductance  : thermal conductance of 1 and 2 stage of TEC
%     .Voltage: input electrical voltage [V]
%     .Current: input electrical current [A]
%  opt   - (i integer scalar) optional running mode
%                             = 0 (default) 采用参考文献[1]的方法计算TEC性能参数
%                             = 1 采用参考文献[2]的方法计算TEC性能参数
%
%   by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
%  References
%  [1] Xuan X C, et al. Cryogenics, 2002, 42: 273-278.
%  [2] Huang B J, et al. International Journal of Refrigeration 2000, 23(3): 208-218.
%
%%  function body
% intialize
switch nargin
    case(2)
        opt = 0;
    case(3)
        if opt ~= 0 && opt ~= 1
            prompt = sprintf('Unknown specified running mode of %d for TE_MaterialProp()', opt);
            TE_log(prompt, 1);
            return
        end
end
results = zeros(0,3);
% parameters of property calculation,
switch opt
    case(0)
        % ref. to eqs.(8)-(10) in [1]
        params = [2.2224e-5,  9.306e-7, -9.905e-10; ...
                   5.112e-7,  1.634e-8,  6.279e-11; ...
                   6.2605  , -2.777e-2,  4.131e-5 ];
        % temperature array
        Ts = [1, T_avg, T_avg^2];               
    case(1)
        % 获取计算a R K的参数
        params = TEC.Parameters;
        TEC.NumRatio = 0; % 将r值设定为0，即按单层TEC结构计算
        % 按计算参数设定温度向量
        Ts = zeros(1, size(params, 2));
        for i = 1:size(params, 2)
            Ts(i) = 1*T_avg^(i-1);
        end
end

% calculate the properties of thermocouple
for i = 1:size(params, 1)
    results(i) = dot(params(i,:), Ts); % results = [alpha rho kappa]
end
%% 输出结果
switch opt
    case(0)
        TEC.SeebeckCoefficient = (results(1)+results(1));
        TEC.ElecConductance    = (results(2)+results(2))/TEC.GeomFactor;
        TEC.ThermConductance   = (results(3)+results(3))*TEC.GeomFactor;
    case(1)
        TEC.SeebeckCoefficient = results(1);
        TEC.ElecConductance    = results(2);
        TEC.ThermConductance   = results(3);
end
%
end

