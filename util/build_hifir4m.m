function build_hifir4m(force)
% script for building HIFIR4M

if nargin < 1; force = false; end
mods = {'mex/hifir4m_mex', ...
    'matlab/private/hifir4m_ijv2crs', ...
    'matlab/private/hifir4m_isint64'};
    
if isoctave
    mexCmd = 'mmex';
    sysLibs = ' -llapack -lblas';
else
    mexCmd = 'mex';
    if ispc
        sysLibs = ' libmwlapack.lib libmwblas.lib';
    else
        sysLibs = ' -lmwlapack -lmwblas';
    end
end

for m = 1:length(mods)
    md = mods{m};
    src = relativepath(fullfile(hifir4m_root, [md '.cpp']));
    mx = relativepath(fullfile(hifir4m_root, [md '.' mexext]));
    if ~force && exist(mx, 'file') && ~isnewer(src, mx); continue; end
    % assume GCC openmp
    cmd = [mexCmd ' ' ...
        'LDFLAGS="$LDFLAGS -fopenmp" ' ... % OpenMP linker flag
        'CXXFLAGS="$CXXFLAGS -m64 -march=native -O3 -std=c++11 ' ...
        '-ffast-math -fcx-limited-range -fopenmp" ' ... % C++11/OpenMP compiler
        '-I' relativepath(fullfile(hifir4m_root, 'hifir', 'src')) ' ' ... % include
        '-O -output ' mx ' ' src sysLibs];  % link to system libraries
    disp(cmd);
    eval(cmd);
end
end

function flag = isnewer(f1, f2)
    d1 = dir(f1);
    d2 = dir(f2);
    flag = datenum(d1.date) >= datenum(d2.date);
end
