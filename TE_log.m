%% 日志记录
%  功能
%  1. 命令行输出当前日志
%  2. 用全局表变量TE_LogData存储日志记录
%
% by Dr. Guan Guoqiang @ SCUT on 2020-04-05
%
function TE_log(content, status, opt)
% 函数调用输入、输出参数说明
% content - (i string) 单行日志内容字段
% status  - (i integer scalar) 可选当前运行状态代码
%                             0 - （缺省值）[INFO]
%                             1 - [ERROR]
%                             2 - [WARNING]
% opt     - (i integer scalar) 可选程序操作代码
%                             0 - (缺省值) 显示当前记录并添加到日志
%                             1 - 添加当前记录到日志变量但不显示

%% 初始化
switch nargin
    case(1)
        status = 0;
        opt = 0;
    case(2)
        opt = 0;
    case(3)
        if opt ~= 0 && opt ~= 1
            TE_log('Unknown input argument in TE_log()', 1);
            return
        end
end
global TE_LogData

%% 生成日志记录
% 日志时间
log_str.datetime = datetime(now, 'ConvertFrom', 'datenum');
% 日志内容
log_str.content = content;
% 操作状态
switch status
    case(0)
        log_str.status = '[INFO]';
    case(1)
        log_str.status = '[ERROR]';
    case(2)
        log_str.status = '[WARNING]';
    otherwise
        log_str.status = '[UNKNOWN]';
end
output = struct2table(log_str, 'AsArray', 1);

%% 输出
switch opt
    case(0)
        % 命令行显示当前日志记录
        fprintf('%s | %s %s\n', datestr(output.datetime, 31), output.status{:}, output.content{:})
    case(1)
end
% 转换日志记录为表并添加到全局表变量TE_LogData
if exist('TE_LogData', 'var') == 1
    TE_LogData = [TE_LogData; output];
else
    TE_LogData = output;
end

end