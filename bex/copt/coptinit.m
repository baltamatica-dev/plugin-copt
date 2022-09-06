function ret = coptinit()
%coptinit - 初始化 copt 插件
%
% Usage:
%   ret = coptinit();
%
% 执行此脚本以初始化 copt 插件.
% 返回值
%   - 为 0 说明初始化成功
%   - 其他值说明初始化失败
%
% 初始化函数可以多次执行.

PLUGIN_NAME = "copt";
ret = -1;
unload_plugin(PLUGIN_NAME);
ret = load_plugin(PLUGIN_NAME);

% 加载当前路径下脚本
addpath(".");
addpath("export");


% 判断插件是否加载成功
if (0==ret)
    disp(PLUGIN_NAME + " 插件加载成功.");
else
    disp(PLUGIN_NAME + " 插件失败.");
    disp("请参考以上出现的报错信息进行修正.");
end % if (0==ret)

end % SPDX-License-Identifier: MIT