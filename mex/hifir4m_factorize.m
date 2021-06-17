function varargout = hifir4m_factorize(dbase, A, varargin)
%HIFIR4M_FACTORIZE - Factorize HIFIR preconditioner
%
% Syntax:
%   hifir4m_factorize(dbase, A)
%   hifir4m_factorize(dbase, A, opts)
%   hifir4m_factorize(dbase, A, fieldName, fieldValue, ...)
%   t = hifir4m_factorize(___)
%   [t, fac_info] = hifir4m_factorize(___)
%
% Description:
%   HIFIR4M_FACTORIZE is the driver call for factorizing a squared sparse
%   matrix with HIFIR. It has the above different syntaxes, which will be
%   addresses shortly.
%
%   hifir4m_factorize(dbase, A) simply factorizes a given matrix with a
%   pre-constructed internal database "dbase" and a squared sparse matrix A.
%   This syntax uses the default options.
%
%   hifir4m_factorize(dbase, A, opts) is like above but with a customized
%   option structured.
%
%   hifir4m_factorize(dbase, A, fieldName, fieldValue, ...) allows one to
%   implicitly pass in parameters without explicitly dealing with an option
%   structure.
%
%   [t, fac_info] = hifir4m_factorize(___) indicates that the function has
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
%   For those who are directly working with CRS in MATLAB, they can wrap their
%   CRS matrix into a structure with fields:
%       row_ptr - starting position array (offset)
%       col_ind - column index array
%       val     - value array
%   Notice that for structure input, the row_ptr and col_ind must be int32
%
% Examples:
%   To simply factorize a sparse matrix, we can
%       >> dbase = hifir4m_initialize;
%       >> A = sprand(10, 10, 0.5);
%       >> hifir4m_factorize(dbase, A);
%
%   If you are interested in the timing
%       >> t = hifir4m_factorize(dbase, A);
%       >> disp(t);
%
%   Supply your own control options
%       >> hifir4m_factorize(dbase, A, 'symm_pre_lvls', 2); % 2 level symm
%
%   Pass the matrix via MATFILE
%       >> A = sprand(10, 10, 0.5);
%       >> save test.mat A
%       >> clear A
%       >> hifir4m_factorize(dbase, {'test.mat', 'A'});
%
%   Using struct
%       >> A.row_ptr = int32(...);
%       >> A.col_ind = int32(...);
%       >> A.val = ...;
%       >> hifir4m_factorize(dbase, A);
%
% See Also:
%   HIFIR4M_INITIALIZE, HIFIR4M_CREATE_PARAMS, HIFIR4M_FGMRES

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: AGPLv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

if iscell(A)
    assert(length(A) == 2);
    t = tic; A = getfield(load(A{1}, A{2}), A{2}); t = toc(t);
    fprintf(1, 'HIFIR4M factorization I/O time is %.4gs\n', t);
end
if ~isempty(varargin) && isstruct(varargin{1})
    opts = varargin{1};  % ignore whatever goes after the first one
else
    opts = hifir4m_create_params(varargin{:});
end
if nargin < 3 || isempty(opts);  end
% convert A to zero-based CRS
if issparse(A); A = hifir4m_sp2crs(A); end
A = hifir4m_ensure_int(A);
[varargout{1:nargout}] = hifir4m_mex(HIFIR4M_FACTORIZE, dbase, ...
    A.row_ptr, A.col_ind, A.val, opts);

%-------------------------- END MAIN CODE -------------------------------%
end
