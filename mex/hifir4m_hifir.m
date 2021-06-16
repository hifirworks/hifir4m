function varargout = hifir4m_hifir(dbase, varargin)
%HIFIR4M_HIFIR - Accessing inv(M) using iterative refinement
%
% Syntax:
%   x = hifir4m_hifir(dbase, A, b, N)
%   x = hifir4m_hifir(dbase, A, b, N, trans)
%   x = hifir4m_hifir(dbase, A, b, N, trans, r)
%   [x, t] = hifir4m_hifir(___)
%
% Description:
%   HIFIR4M_HIFIR allows one to directly access the preconditioner in a
%   iterative-refinement fashion.
%
%   x = hifir4m_hifir(dbase, A, b, N) iterative refinements with M^g and A
%   with N iterations.
%
%   x = hifir4m_hifir(dbase, A, b, N, trans, r) iteratively refines M^g and A
%   with operation trans (i.e., true for A' and (M^g)') and final Schur
%   complement dimension r.
%
%   [x, t] = hifir4m_hifir(___) can potentially report the
%   overhead-free wall-clock time.
%
% See Also:
%   HIFIR4M_SOLVE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: AGPLv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

assert(length(varargin) >= 3, 'inputs must be at least A, b, and N');
b = varargin{2};
N = varargin{3}; N = N(1);
% transpose/Hermitian flag
trans = false;
if length(varargin) > 3; trans = logical(varargin{4}); trans = trans(1); end
% final Schur complement dimension
r = -1;
if length(varargin) > 4; r = varargin{5}; r = r(1); end
A = varargin{1};
if issparse(A); A = hifir4m_sp2crs(A); end
A = hifir4m_ensure_int(A);
[varargout{1:nargout}] = hifir4m_mex(HIFIR4M_M_SOLVE, dbase, b, ...
    A.row_ptr, A.col_ind, A.val, N, trans, r);

%-------------------------- END MAIN CODE -------------------------------%
end