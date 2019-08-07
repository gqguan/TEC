%% investigate the COP vs. Qc/N/(r+1) for two-stage serial TEC
% 
% by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
clear;
% set the waitbar
ui_waitbar = waitbar(0, 'Setting ... ');
% set parameters
N = 127;
Th = 300; 
Tc = 250;
r = 5;
% 计算电流上下边界
IBound = TE_Current(Th, Tc, N, r, 1);
IMax = max(IBound);
IMin = min(IBound);
dI = (IMax-IMin)/100;
% initialize
I = IMin:dI:IMax;
COP = zeros(size(I));
XCoordinate = zeros(size(I));
%
waitbar(0, ui_waitbar, 'Calculating ... ')
% fprintf('Calculating ... [%d]: ', length(I));
for i = 1:length(I)
    wb_str = sprintf('Calculating ... %4.1f %%', (i/length(I)*100));
%     wb_str = ['Calculating ...', num2str(i/length(I)*100),'%'];
    waitbar(i/length(I), ui_waitbar, wb_str)    
    [Q, ~, ~, ~] = TE_Heat(Th, Tc, I(i), N, r);
    COP(i) = Q(2)/(Q(1)-Q(2));
    XCoordinate(i) = Q(2)/N/(r+1);
end
delete(ui_waitbar);
fprintf(' completed \n');