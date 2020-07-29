function hilu = hilucsi4m_export(dbase, varargin)
%HILUCSI4M_EXPORT - Export internal data to MATLAB
%
% Syntax:
%   hilu = hilucsi4m_export(dbase)
%   hilu = hilucsi4m_export(dbase, Options)
%
% Description:
%   HILUCSI4M_EXPORT exports the internal data inside C++ to MATLAB
%
%   hilu = hilucsi4m_export(dbase) exports data in all levels to MATLAB,
%   where hilu is a cell array, whose length is the number of levels.
%   Each sparse level will contain a struct of "L_B", "D_B, "U_B", "E",
%   "F", "s_row", "s_col", "p", "p_inv", "q", and "q_inv" total 11
%   attributes. For the last dense level, it will output the dense matrix,
%   and the user need to factorize it manullly by either LU or QR.
%
%   hilu = hilucsi4m_export(dbase, Options) with customized options:
%       1. 'format': {'sparse', 'economic'}, the default is `sparse`. For
%           'economic', it will output struct of CRS_MATRIX.
%       2. 'destroy': boolean flag, if true, then it will destroy the
%           internal data level by level while exporting it to MATLAB.
%           The default value is false. (Advanced usage!)
%
% Examples:
%   To get all levels in sparse format
%       >> hilu = hilucsi4m_export(dbase);
%       >> for lvl = 1:size(hilu,1); disp(hilu{i}); end
%
%   To get with economic format, do the following
%       >> hilu = hilucsi4m_export(dbase, 'format', 'economic');
%
% See Also:
%   HILUCSI4M_FACTORIZE, LU, QR

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: GLPv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

p = inputParser;
addOptional(p, 'format', 'sparse', @(x) ismember(x, {'sparse', 'economic'}));
addOptional(p, 'destroy', false, @(x) isscalar(x));
parse(p, varargin{:});  % parse
opts = p.Results;

hilu = hilucsi4m_mex(HILUCSI4M_EXPORT_DATA, dbase, opts.destroy);

if strcmp(opts.format, 'sparse')
    if isempty(which('crs_2sparse'))
        fprintf(2, 'converting to sparse requires crs_2sparse from numgeom\n');
        return;
    end
    for i = 1:size(hilu,1)
        if isstruct(hilu{i})
            % sparse level
            n = double(hilu{i}.L_B.nrows);
            hilu{i}.L_B = crs_2sparse(hilu{i}.L_B);
            hilu{i}.L_B = hilu{i}.L_B + speye(n);
            hilu{i}.D_B = spdiags(hilu{i}.D_B, 1, n, n);
            hilu{i}.U_B = crs_2sparse(hilu{i}.U_B);
            hilu{i}.U_B = hilu{i}.U_B + speye(n);
            hilu{i}.E = crs_2sparse(hilu{i}.E);
            hilu{i}.F = crs_2sparse(hilu{i}.F);
        end
    end
end

%-------------------------- END MAIN CODE -------------------------------%
end