function startup
% Startup script for HILUCSI4M

addpath(genpath(hilucsi4m_root));
rmpath(genpath([hilucsi4m_root '/hilucsi']));
rmpath(genpath([hilucsi4m_root '/examples']));
if exist(fullfile(hilucsi4m_root, '.git'), 'dir')
    rmpath(genpath(fullfile(hilucsi4m_root, '.git')));
end
if exist(fullfile(hilucsi4m_root, '.vscode'), 'dir')
    rmpath(genpath(fullfile(hilucsi4m_root, '.vscode')));
end

% build mex
build_hilucsi4m;
end