function [Q, a, R, K] = TE_Heat(Th, Tc, I, N, r)
%% calculate the heats of TEC
%  notes of I/O arguments
%  Th - (i double scalar) hot-side temperature [K]
%  Tc - (i double scalar) cold-side temperature [K]
%  I  - (i double scalar) serial electric current flowing through the two
%                         stages [A]
%  N  - (i double scalar) number of thermocouples in the first stage
%  r  - (i double scalar) number ratio of thermocouples in two stages
%  Q  - (o double array(2)) heats flowing out/in the hot/cold side of TEC
%  a  - (o double array(2)) Seebeck coefficient of 1 and 2 stage of TEC
%  R  - (o double array(2)) electrical resistance of 1 and 2 stage of TEC
%  K  - (o double array(2)) thermal conductance of 1 and 2 stage of TEC
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
syms Qh Qc Tm;
% initialize
dTm = 1e5;
direction_Tm = -1;
% % set values
% Th = 320;
% Tc = 290;
% I  = 0.88;
% N  = 127;
% r  = 1;
% 试差法计算Tm
Tm_val = (Th+Tc)/2; % Tm初值
while dTm > 1e-5
    dTm0 = dTm;
    % calculate a1 R1 K1
    [a1, R1, K1] = TE_MaterialProp((Th+Tm_val)/2, 0.0015);
    % calculate a2 R2 K2
    [a2, R2, K2] = TE_MaterialProp((Tc+Tm_val)/2, 0.0015);
    % define heat balance equations
    eq(2) = Qh == (I*a1*Th+I^2*R1/2-K1*(Th-Tm))*N*r;
    eq(3) = (I*a1*Tm-I^2*R1/2-K1*(Th-Tm))*r == I*a2*Tm-I^2*R2/2-K2*(Tm-Tc);
    eq(4) = Qc == (I*a2*Tc-I^2*R2/2-K2*(Tm-Tc))*N;
    % get the junction temperature between stages Tm
    dTm = abs(eval(solve(eq(3), Tm))-Tm_val);
    if (dTm > dTm0)
        direction_Tm = -direction_Tm;
    end
    Tm_val = Tm_val+direction_Tm*dTm;
end
% 吸、放热量
Q(1) = eval(solve(subs(eq(2), Tm, Tm_val), Qh));
Q(2) = eval(solve(subs(eq(4), Tm, Tm_val), Qc));
% 输出热电偶参数
a = [a1, a2];
R = [R1, R2];
K = [K1, K2];
%
end
