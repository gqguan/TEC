%% 公式推导
%
% by Dr. Guan Guoqiang @ SCUT on 2022-05-07

%% 初始化
clear

%% 推导部分回流、稳态操作时DCMD系统料液侧吸热量
syms WF TF0 Q1 Q2 % 待求变量
syms W1 TF1 R % 操作参数
syms WP QM TMF TF2 TMP TP2 % 模拟求得参数
syms cp1 cp2 % 物性参数
syms T0 % 环境参数
% 对于笔记2022/5/6所示系统（classical、extTEHP和permTEHP）
% 对于回流混合单元CV1（见笔记2022/5/16），能量平衡方程
WR1 = R*(WF-WP); % 回流比
cv1EEq = WF*cp1*T0 + WR1*cp1*TF2 == W1*cp1*TF0; % classical或extTEHP
cv1aEEq = WF*cp1*T0 + WR1*cp1*TF2 == W1*cp1*TF1; % feedTEHP或permTEHP1
cv1bEEq = WF*cp1*T0 + WR1*cp1*TF0 == W1*cp1*TF1; % permTEHP2
% 对于DCMD料液侧CV2（见笔记2022/5/16），能量平衡方程（即为系统内每个单元都进行能量平衡综合的结果）
cv2EEq = WF*cp1*T0 + Q1 == (WF-WP)*cp1*TF2+WP*cp1*TMF+QM; % classical、extTEHP、feedTEHP或permTEHP1/2
% 对于料液加热单元（见笔记2022/5/16），能量平衡方程
cv3EEq = W1*cp1*TF0 + Q1 == W1*cp1*TF1; % classical、extTEHP、feedTEHP或permTEHP1
cv3bEEq = WR1*cp1*TF0 == WR1*cp1*TF2+Q1; % permTEHP2
% 对于DCMD渗透侧CV4（见笔记2022/5/16）
cv4EEq = Q2+WP*cp2*TP2 == QM+WP*cp2*TMP;
% 对于DCMD系统CV5（见笔记2022/5/18）
cv5EEq = WF*cp1*T0+(Q1-Q2) == WP*cp2*TP2+(WF-WP)*cp1*TF2;

% 联立求解
sol1 = solve([cv1EEq,cv2EEq,cv3EEq],[WF,TF0,Q1]);
disp('对于classical、extTEHP和permTEHP1')
disp('WF = ')
disp(sol1.WF)
disp('TF0 = ')
disp(sol1.TF0)
disp('Q1 = ')
disp(sol1.Q1)

sol1a = solve([cv1aEEq,cv2EEq],[WF,Q1]);
disp('对于feedTEHP')
disp('WF = ')
disp(sol1a.WF)
disp('Q1 = ')
disp(sol1a.Q1)

sol1b = solve([cv1bEEq,cv2EEq,cv3bEEq],[WF,TF0,Q1]);
disp('对于permTEHP2')
disp('WF = ')
disp(sol1b.WF)
disp('TF0 = ')
disp(sol1b.TF0)
disp('Q1 = ')
disp(sol1b.Q1)

%% 稳态操作时给定DCMD系统料液侧吸热量计算回流比及料液处理量
% 联立求解
sol2 = solve([cv1EEq,cv2EEq,cv3EEq],[WF,R,TF0]);
% 显示结果
disp('classical或extTEHP：R = ')
disp(sol2.R)
