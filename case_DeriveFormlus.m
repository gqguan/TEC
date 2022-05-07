%% 推导部分回流、稳态操作时DCMD系统料液侧吸热量
%
% by Dr. Guan Guoqiang @ SCUT on 2022-05-07

%% 初始化
clear
syms WF TF0 Q1 % 待求变量
syms W1 TF1 R % 操作参数
syms WP QM TMF TF2 % 模拟求得参数
syms cp % 物性参数
syms T0 % 环境参数

%% 对于笔记2022/5/6所示系统
% 对于回流混合单元CV1，能量平衡方程
WR1 = R*(WF-WP);
cv1EEq = WF*cp*T0 + WR1*cp*TF2 == W1*cp*TF0;
% 对于DCMD料液侧CV2，能量平衡方程（即为系统内每个单元都进行能量平衡综合的结果）
cv2EEq = WF*cp*T0 + Q1 == (WF-WP)*cp*TF2+WP*cp*TMF+QM;
% 对于料液加热单元，能量平衡方程
cv3EEq = W1*cp*TF0 + Q1 == W1*cp*TF1;
% 对于膜组件料液侧，能量平衡方程
cv4EEq = W1*cp*TF1 == (WR1+(WF-WP))*cp*TF1+WP*cp*TMF+QM;
% 联立求解
solEEqSys = solve([cv1EEq,cv3EEq,cv4EEq],[WF,TF0,Q1]);

disp(solEEqSys.Q1)

