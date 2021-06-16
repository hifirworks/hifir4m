function hilu = hifir4m_export(dbase, varargin)
%HIFIR4M_EXPORT - Export internal data to MATLAB
%
% Syntax:
%   hilu = hifir4m_export(dbase)
%   hilu = hifir4m_export(dbase, Options)
%
% Description:
%   HIFIR4M_EXPORT exports the internal data inside C++ to MATLAB
%
%   hilu = hifir4m_export(dbase) exports data in all levels to MATLAB,
%   where hilu is a cell array, whose length is the number of levels.
%   Each sparse level will contain a struct of "L_B", "D_B, "U_B", "E",
%   "F", "s_row", "s_col", "p", "p_inv", "q", and "q_inv" total 11
%   attributes. For the last dense level, it will output the dense matrix,
%   and the user need to factorize it manullly by either LU or QR.
%
%   hilu = hifir4m_export(dbase, Options) with customized options:
%       1. 'format': {'sparse', 'economic'}, the default is `sparse`. For
%           'economic', it will output struct of CRS_MATRIX.
%       2. 'destroy': boolean flag, if true, then it will destroy the
%           internal data level by level while exporting it to MATLAB.
%           The default value is false. (Advanced usage!)
%
% Examples:
%   To get all levels in sparse format
%       >> hilu = hifir4m_export(dbase);
%       >> for lvl = 1:size(hilu,1); disp(hilu{i}); end
%
%   To get with economic format, do the following
%       >> hilu = hifir4m_export(dbase, 'format', 'economic');
%
% See Also:
%   HIFIR4M_FACTORIZE, LU, QR

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: AGPLv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

p = inputParser;
addOptional(p, 'format', 'sparse', @(x) ismember(x, {'sparse', 'economic'}));
addOptional(p, 'destroy', false, @(x) isscalar(x));
parse(p, varargin{:});  % parse
opts = p.Results;

hilu = hifir4m_mex(HIFIR4M_EXPORT_DATA, dbase, opts.destroy);

if strcmp(opts.format, 'sparse')
    for i = 1:size(hilu,1)
        if isstruct(hilu{i})
            % sparse level
            n = double(hilu{i}.L_B.nrows);
            hilu{i}.L_B = crsSparse(hilu{i}.L_B);
            hilu{i}.L_B = hilu{i}.L_B + speye(n);
            hilu{i}.D_B = spdiags(hilu{i}.D_B, 1, n, n);
            hilu{i}.U_B = crsSparse(hilu{i}.U_B);
            hilu{i}.U_B = hilu{i}.U_B + speye(n);
            hilu{i}.E = crsSparse(hilu{i}.E);
            hilu{i}.F = crsSparse(hilu{i}.F);
        end
    end
end

%-------------------------- END MAIN CODE -------------------------------%
end

function row_ind = crsRowind(row_ptr, col_ind)
row_ind = zeros(size(col_ind),class(col_ind));

nrows = int32(length(row_ptr))-1;
for i=1:nrows
    for j = row_ptr(i) : row_ptr(i+1) - 1
        row_ind(j) = i;
    end
end
end

function B = crsSparse(A)
row_ind = crsRowind(A.row_ptr, A.col_ind);
B = sparse(double(row_ind), double(A.col_ind), A.val, double(A.nrows), ...
    double(A.ncols));
end