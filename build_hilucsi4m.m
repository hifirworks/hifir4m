function build_hilucsi4m(force)
% script for building hilucsi4m

if nargin < 1; force = false; end
mods = {'mex/hilucsi4m_mex', ...
    'matlab/private/hilucsi4m_ijv2crs'};
%{
if system('git submodule update --init hilucsi')
    fprintf(2, 'Warning! Failed to update submodule HILUCSI\n');
end
%}
for m = 1:length(mods)
    md = mods{m};
    src = fullfile(hilucsi4m_root, [md '.cpp']);
    mx = fullfile(hilucsi4m_root, [md '.' mexext]);
    if ~force && exist(mx, 'file') && ~isnewer(src, mx); continue; end
    % assume GCC openmp
    cmd = ['mex ' ...
        'LDFLAGS="$LDFLAGS -fopenmp" ' ... % OpenMP linker flag
        'CXXFLAGS="$CXXFLAGS -m64 -march=native -O3 -std=c++11 -fopenmp" ' ... % C++11/OpenMP compiler
        '-I' fullfile(hilucsi4m_root, 'hilucsi', 'src') ' ' ... % include
        '-v -O -R2018a -output ' mx ' ' src ...
        ' -lmwblas -lmwlapack' ...  % link to MATLAB Lapack/BLAS
        ' -lmwservices'];  % link to internal ioFlush
    try
        disp(cmd);
        eval(cmd);
    catch
        error('Error during compilcation with err %s.', lasterr);
    end
end
end

function flag = isnewer(f1, f2)
    d1 = dir(f1);
    d2 = dir(f2);
    flag = datenum(d1.date) >= datenum(d2.date);
end
