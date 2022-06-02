%% calculate the heats of TEC
%  notes of I/O arguments
%  Th  - (i double scalar) hot-side temperature [K]
%  Tc  - (i double scalar) cold-side temperature [K]
%  I   - (i double scalar) serial electric current flowing through the two
%                         stages [A]
%  TEC - (i struc) struc variable
%     .NumTC     : Number of thermocouples in TEC
%     .NumRatio  : ratio of thermocouples in the 1-stage TEC to those in
%                    the 2-stage TEC
%     .GeomFactor: geometry factor of thermcouples in TEC [m]
%     .SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC
%     .ElecConductance   : electrical conductance of 1 and 2 stage of TEC
%     .ThermConductance  : thermal conductance of 1 and 2 stage of TEC
%     .Voltage: input electrical voltage [V]
%     .Current: input electrical current [A]
%  opt1 - (i optional scalar) calculation method
%        = 0 (default) 按参考文献[1]计算吸放热量
%        = 1           按参考文献[2]计算吸放热量
%  opt2 - (i optional scalar) given electrical current or voltage
%        = 0 (default) 给定电流计算吸放热量
%        = 1           给定电压计算吸放热量
%  Q   - (o double array(2)) heats flowing out/in the hot/cold side of TEC
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
%  2019-08-07: add case 0 for calculating the heats in one-stage TEC
%  2019-08-10: update according to the new TE_Tm()
%  2019-08-19: change struct(TEC) to include the U and I
%  2020-04-03: 添加运行模式设定
%
%  References
%  [1] Xuan X C, et al. Cryogenics, 2002, 42: 273-278.
%  [2] Huang B J, et al. International Journal of Refrigeration 2000, 23(3): 208-218.
%  [3] Chen L et al. Applied Energy, 2008, 85: 641-649 [doi:10.1016/j.apenergy.2007.10.005]
%
function [Q, TEC] = TE_Heat(Th, Tc, TEC, opt1, opt2)
%% 初始化
% 输入参数检查
if ~exist('opt1','var')
    opt1 = 0; % opt1缺省值
end
if ~exist('opt2','var')
    opt2 = 0; % opt2缺省值
end
% 指定运行模式
if opt1 ~= 0 && opt1 ~= 1
    prompt = sprintf('Unknown specified running mode of %d for TE_Heat() input argument opt1', opt1);
    TE_log(prompt, 1);
    return
end
% Calculate parameters of thermocouples
TEC = TE_MaterialProp((Th+Tc)/2,TEC,opt1);
a = TEC.SeebeckCoefficient;
R = TEC.ElecConductance;
K = TEC.ThermConductance;
if all([a R K] == 0)
    Q = [0 0];
    return
end
switch opt2
    case(0)
        I = TEC.Current;
    case(1)
        U = TEC.Voltage;
        I = (U-TEC.SeebeckCoefficient*(Th-Tc))/TEC.ElecConductance; % eq.(1) ref. [2]
    otherwise
        prompt = sprintf('Unknown specified running mode of %d for TE_Heat() input argument opt2', opt2);
        TE_log(prompt, 1);
        return        
end
%%
switch opt1
    case(0)     
        % Calculate the absorbed and released heats
        switch TEC.NumRatio
            case 0 % 单层结构
                N0 = TEC.NumTC;
                Q(1) = (I*a*Th+I^2*R/2-K*(Th-Tc))*N0;
                Q(2) = (I*a*Tc-I^2*R/2-K*(Th-Tc))*N0;          
%                 % 输出TEC性能参数
%                 TEC.SeebeckCoefficient = a;
%                 TEC.ElecConductance = R;
%                 TEC.ThermConductance = K;                
            otherwise % 两层结构
                m = round(TEC.NumTC*(TEC.NumRatio/(1+TEC.NumRatio)));
                n = TEC.NumTC-m;
                k1 = TEC.HTCoefficient;
                k2 = TEC.HTCoefficient;
                F1 = TEC.HTArea/(1+TEC.NumRatio);
                F2 = TEC.HTArea/(1+TEC.NumRatio)*TEC.NumRatio;
                % 吸、放热量（参考TE_Sandbox1.m中eq.12和13的推导结果）
                T1 = (2*I^2*K^2*R*m*n^2 + 2*I^2*K^2*R*m^2*n + I^4*R*a^2*m*n^2 - I^4*R*a^2*m^2*n + 2*F2*I^2*K*R*k2*n^2 + F2*I^3*R*a*k2*n^2 + 3*I^3*K*R*a*m*n^2 + I^3*K*R*a*m^2*n - 2*F1*I^2*Th*a^2*k1*m^2 + 2*F2*K^2*Tc*k2*m*n + 2*F1*K^2*Th*k1*m*n + 2*F2*I^2*K*R*k2*m*n - F2*I^3*R*a*k2*m*n + 2*F1*I^2*Th*a^2*k1*m*n + 2*F1*F2*K*Th*k1*k2*m + 2*F1*F2*K*Th*k1*k2*n - 2*F1*F2*I*Th*a*k1*k2*m + 2*F1*F2*I*Th*a*k1*k2*n + 4*F1*I*K*Th*a*k1*m*n)/(2*(I^3*a^3*m^2*n - I^3*a^3*m*n^2 - F1*I^2*a^2*k1*m^2 - F2*I^2*a^2*k2*n^2 - I^2*K*a^2*m*n^2 - I^2*K*a^2*m^2*n + F1*K^2*k1*m*n + F2*K^2*k2*m*n + F1*I^2*a^2*k1*m*n + F2*I^2*a^2*k2*m*n + F1*F2*K*k1*k2*m + F1*F2*K*k1*k2*n - F1*F2*I*a*k1*k2*m + F1*F2*I*a*k1*k2*n + 2*F1*I*K*a*k1*m*n - 2*F2*I*K*a*k2*m*n));
                T2 = (2*I^2*K^2*R*m*n^2 + 2*I^2*K^2*R*m^2*n - I^4*R*a^2*m*n^2 + I^4*R*a^2*m^2*n + 2*F1*I^2*K*R*k1*m^2 - F1*I^3*R*a*k1*m^2 - I^3*K*R*a*m*n^2 - 3*I^3*K*R*a*m^2*n - 2*F2*I^2*Tc*a^2*k2*n^2 + 2*F2*K^2*Tc*k2*m*n + 2*F1*K^2*Th*k1*m*n + 2*F1*I^2*K*R*k1*m*n + F1*I^3*R*a*k1*m*n + 2*F2*I^2*Tc*a^2*k2*m*n + 2*F1*F2*K*Tc*k1*k2*m + 2*F1*F2*K*Tc*k1*k2*n - 2*F1*F2*I*Tc*a*k1*k2*m + 2*F1*F2*I*Tc*a*k1*k2*n - 4*F2*I*K*Tc*a*k2*m*n)/(2*(I^3*a^3*m^2*n - I^3*a^3*m*n^2 - F1*I^2*a^2*k1*m^2 - F2*I^2*a^2*k2*n^2 - I^2*K*a^2*m*n^2 - I^2*K*a^2*m^2*n + F1*K^2*k1*m*n + F2*K^2*k2*m*n + F1*I^2*a^2*k1*m*n + F2*I^2*a^2*k2*m*n + F1*F2*K*k1*k2*m + F1*F2*K*k1*k2*n - F1*F2*I*a*k1*k2*m + F1*F2*I*a*k1*k2*n + 2*F1*I*K*a*k1*m*n - 2*F2*I*K*a*k2*m*n));
                Q(1) = F1*k1*(T1-Th);
                Q(2) = F2*k2*(Tc-T2);
%                 Tm = ((m+n)*I^2*R/2+n*K*Th+m*K*Tc)/(a*I*(n-m)+K*(m+n));
%                 Q(1) = (I*a*Th+I^2*R/2-K*(Th-Tm))*n;
        end
     case(1)
%         % TEC性能参数
%         TEC = TE_MaterialProp((Th+Tc)/2, TEC, opt1);
%         a = TEC.SeebeckCoefficient;
%         R = TEC.ElecConductance;
%         K = TEC.ThermConductance;
        % 计算吸放热量
        Q(1) = (I*a*Th+I^2*R/2-K*(Th-Tc)); % 单位：W（下同）
        Q(2) = (I*a*Tc-I^2*R/2-K*(Th-Tc));
end
%
end
