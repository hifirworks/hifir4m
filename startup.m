function startup
% Startup script for HILUCSI4M

addpath(genpath(hilucsi4m_root));
rmpath(genpath([hilucsi4m_root '/hilucsi']));
rmpath(genpath([hilucsi4m_root '/examples']));

% build mex
build_hilucsi4m;
end