function build_hifir4m(varargin)
% script for building hifir4m for MATLAB and GNU Octave

if nargin < 1
    force = false;
end

mods = {'mex/hifir4m_mex', ...
    'mex/private/hifir4m_ijv2crs', ...
    'mex/hifir4m_isint64', ...
    'mex/hifir4m_version'};

if isoctave
    mexCmd = 'mmex -O';
    sysFlags = ' -llapack -lblas';
else
    mexCmd = 'mex -O';
    if ispc
        sysFlags = ' -DHIF_FC=1 -R2018a LINKLIBS=''-llibmwmathutil $LINKLIBS''';
    else
        sysFlags = ' -R2018a -lmwlapack -lmwblas -lmwservices';
    end
end

% download hifir
hifir_verion = '0.1.0';
downloaded = false;

% check HIFIR C++ root
if exist(fullfile(hifir4m_root, ['hifir-' hifir_verion]), 'dir')
    path_to_hifir = fullfile(hifir4m_root, ['hifir-' hifir_verion]);
else
    downloaded = true;
    path_to_hifir = fullfile(tempdir, ['hifir-' hifir_verion]);
end

for m = 1:length(mods)
    md = mods{m};
    src = fullfile(hifir4m_root, [md '.cpp']);
    mx = fullfile(hifir4m_root, [md '.' mexext]);
    if ~force && exist(mx, 'file') && ~isnewer(src, mx); continue; end
    if downloaded && ~exist(path_to_hifir, 'dir')
        path_to_hifir = download_hifir(hifir_verion);
    end
    % assume GCC openmp
    cmd = [mexCmd ' ' ...
        'LDFLAGS="$LDFLAGS -fopenmp" ' ... % OpenMP linker flag
        'CXXFLAGS="$CXXFLAGS -m64 -march=native -O3 -std=c++11 ' ...
        '-ffast-math -fcx-limited-range -fopenmp" ' ... % C++11/OpenMP compiler
        '-I''' fullfile(path_to_hifir, 'src') ''' ' ... % include
        '-output ''' mx ''' ''' src '''' sysFlags];  % link to system libraries
    disp(cmd);
    eval(cmd);
end

% Delete downloaded folder if neccessary
if downloaded && exist(path_to_hifir, 'dir')
    rmdir(path_to_hifir, 's');
end

end

function flag = isnewer(f1, f2)
    d1 = dir(f1);
    d2 = dir(f2);
    flag = datenum(d1.date) >= datenum(d2.date);
end
