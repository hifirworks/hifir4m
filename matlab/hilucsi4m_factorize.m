function varargout = hilucsi4m_factorize(dbase, A, varargin)
%HILUCSI4M_FACTORIZE - Factorize HILUCSI preconditioner
%
% Syntax:
%   hilucsi4m_factorize(dbase, A)
%   hilucsi4m_factorize(dbase, A, opts)
%   hilucsi4m_factorize(dbase, A, fieldName, fieldValue, ...)
%   t = hilucsi4m_factorize(___)
%   [t, fac_info] = hilucsi4m_factorize(___)
%
% Description:
%   HILUCSI4M_FACTORIZE is the driver call for factorizing a squared sparse
%   matrix with HILUCSI. It has the above different syntaxes, which will be
%   addresses shortly.
%
%   hilucsi4m_factorize(dbase, A) simply factorizes a given matrix with a
%   pre-constructed internal database "dbase" and a squared sparse matrix A.
%   This syntax uses the default options.
%
%   hilucsi4m_factorize(dbase, A, opts) is like above but with a customized
%   option structured.
%
%   hilucsi4m_factorize(dbase, A, fieldName, fieldValue, ...) allows one to
%   implicitly pass in parameters without explicitly dealing with an option
%   structure.
%
%   [t, fac_info] = hilucsi4m_factorize(___) indicates that the function has
%   (potentially) two outputs, where the first one *t* is the total factorizing
%   time without MATLAB interpreter overhead. *fac_info* is a structure of
%   some information field of interests resulting from the factorization
%   computation.
%
%   It is worth noting that the variable A can also be a set of length 2,
%   in which the first element is a matfile (by name) containing the target
%   sparse matrix whose variable name is given by the second element in A.
%   The matrix will be cleared once it has been converted into CRS.
%
% Examples:
%   To simply factorize a sparse matrix, we can
%       >> dbase = hilucsi4m_initialize;
%       >> A = sprand(10, 10, 0.5);
%       >> hilucsi4m_factorize(dbase, A);
%
%   If you are interested in the timing
%       >> t = hilucsi4m_factorize(dbase, A);
%       >> disp(t);
%
%   Supply your own control options
%       >> hilucsi4m_factorize(dbase, A, 'symm_pre_lvls', 2); % 2 level symm
%
%   Pass the matrix via MATFILE
%       >> A = sprand(10, 10, 0.5);
%       >> save test.mat A
%       >> clear A
%       >> hilucsi4m_factorize(dbase, {'test.mat', 'A'});
%
% See Also:
%   HILUCSI4M_INITIALIZE, HILUCSI4M_CREATE_OPTIONS, HILUCSI4M_FGMRES

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: GLPv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

via_file = false;
if iscell(A)
    assert(length(A) == 2);
    via_file = true;
    A = getfield(load(A{1}, A{2}), A{2});
end
if ~isempty(varargin) && isstruct(varargin{1})
    opts = varargin{1};  % ignore whatever goes after the first one
else
    opts = hilucsi4m_create_options(varargin{:});
end
if nargin < 3 || isempty(opts);  end
% convert A to zero-based CRS
[rowptr, colind, vals] = hilucsi4m_sp2crs(A);
if via_file; clear A; end
[varargout{1:nargout}] = hilucsi4m_mex(HILUCSI4M_FACTORIZE, dbase, rowptr, ...
    colind, vals, opts);

%-------------------------- END MAIN CODE -------------------------------%
end
