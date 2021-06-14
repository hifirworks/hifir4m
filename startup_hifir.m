function startup_hifir
% Startup script of hifir4m for MATLAB and Octave

addpath(hifir4m_root);
addpath([hifir4m_root '/matlab']);
addpath([hifir4m_root '/util']);

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
