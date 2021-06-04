function startup
% Startup script for HIFIR4M

addpath(genpath(hifir4m_root));
rmpath(genpath([hifir4m_root '/hifir']));
rmpath(genpath([hifir4m_root '/examples']));
if exist(fullfile(hifir4m_root, '.git'), 'dir')
    rmpath(genpath(fullfile(hifir4m_root, '.git')));
end
if exist(fullfile(hifir4m_root, '.vscode'), 'dir')
    rmpath(genpath(fullfile(hifir4m_root, '.vscode')));
end

% build mex
build_hifir4m;
end
