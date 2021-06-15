function startup_hifir
% Startup script of hifir4m for MATLAB and Octave

addpath(hifir4m_root);
addpath([hifir4m_root '/util']);
addpath([hifir4m_root '/api']);
addpath([hifir4m_root '/api/crs']);
addpath([hifir4m_root '/ksp']);
addpath([hifir4m_root '/mex']);

if ~exist('coder.p', 'file')
    addpath([hifir4m_root '/No_coder']);
end

if exist('OCTAVE_VERSION', 'builtin')
    more off;
    warning('off', 'Octave:shadowed-function');
    warning('off', 'Octave:data-file-in-path')
    addpath([hifir4m_root '/util/octave']);
else
    warning('off', 'MATLAB:mex:GccVersion');
    warning('off', 'MATLAB:mex:GccVersion_link');
end

% build mex
build_hifir4m;
end
