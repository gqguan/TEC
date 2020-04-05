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
%  opt - (i optional scalar) running mode
%        = 0 (default) 按参考文献[1]计算吸放热量
%        = 1           按参考文献[2]计算吸放热量
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
%
function [Q, TEC] = TE_Heat(Th, Tc, TEC, opt)
%% 初始化
% 运行模式设定
switch nargin
    case(3)
        % 缺省运行模式
        opt = 0;
    case(4)
        % 指定运行模式
        if opt ~= 0 && opt ~= 1
            prompt = sprintf('Unknown specified running mode of %d for TE_Heat()', opt);
            TE_log(prompt, 1);
            return
        end
end
I = TEC.Current;
%%
switch opt
    case(0)
        % Calculate parameters of thermocouples
        N0 = TEC.NumTC/(TEC.NumRatio+1);      
        % Calculate the absorbed and released heats
        switch TEC.NumRatio
            case 0 % 单层结构
                TEC = TE_MaterialProp((Th+Tc)/2, TEC);
                a = TEC.SeebeckCoefficient;
                R = TEC.ElecConductance;
                K = TEC.ThermConductance;
                Q(1) = (I*a*Th+I^2*R/2-K*(Th-Tc))*N0;
                Q(2) = (I*a*Tc-I^2*R/2-K*(Th-Tc))*N0;
                % 输出TEC性能参数
                TEC.SeebeckCoefficient = a;
                TEC.ElecConductance = R;
                TEC.ThermConductance = K;                
            otherwise % 两层结构
                % Get Tm
                [Tm, TEC] = TE_Tm(Th, Tc, I, TEC);
                a = TEC.SeebeckCoefficient;
                R = TEC.ElecConductance;
                K = TEC.ThermConductance;
                % 吸、放热量
                Q(1) = (I*a(1)*Th+I^2*R(1)/2-K(1)*(Th-Tm))*N0*TEC.NumRatio;
                Q(2) = (I*a(2)*Tc-I^2*R(2)/2-K(2)*(Tm-Tc))*N0;
        end
     case(1)
        % 输入温度应为摄氏度，当输入温度大于200时识别为绝对温度
        if Th > 200
            Th = Th-273.15;
        end
        if Tc > 100
            Tc = Tc-273.15;
        end
        % TEC性能参数
        TEC = TE_MaterialProp((Th+Tc)/2, TEC, opt);
        a = TEC.SeebeckCoefficient;
        R = TEC.ElecConductance;
        K = TEC.ThermConductance;
        % 计算吸放热量
        Q(1) = (I*a*Th+I^2*R/2-K*(Th-Tc)); % 单位：W（下同）
        Q(2) = (I*a*Tc-I^2*R/2-K*(Th-Tc));
end
%
end
