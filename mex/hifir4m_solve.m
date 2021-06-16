function varargout = hifir4m_solve(dbase, varargin)
%HIFIR4M_SOLVER - Accessing inv(M)
%
% Syntax:
%   x = hifir4m_solve(dbase, b)
%   x = hifir4m_solve(dbase, b, trans)
%   x = hifir4m_solve(dbase, b, trans, r)
%   [x, t] = hifir4m_solve(___)
%
% Description:
%   HIFIR4M_SOLVER allows one to directly access the preconditioner,
%   which is useful for developers who want to integrate HIFIR into
%   existing (linear) solver frameworks.
%
%   x = hifir4m_solve(dbase, b) computes M\b given a factorized
%   database.
%
%   x = hifir4m_solve(dbase, b, trans) computes M'\b given a factorized
%   database if trans==true.
%
%   x = hifir4m_solve(dbase, b, trans, r) computes a multilevel triangular
%   solve with operation trans with a customized final Schur complement
%   dimension r.
%
%   [x, t] = hifir4m_solve(___) can potentially report the
%   overhead-free wall-clock time.
%
% See Also:
%   HIFIR4M_FACTORIZE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: AGPLv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

assert(length(varargin) >= 1, 'input(s) must be at least RHS b.');
% transpose/Hermitian flag
trans = false;
if length(varargin) > 1; trans = logical(varargin{2}); trans = trans(1); end
% final Schur complement dimension
r = 0;
if length(varargin) > 2; r = varargin{3}; r = r(1); end
[varargout{1:nargout}] = hifir4m_mex(HIFIR4M_M_SOLVE, dbase, varargin{1}, ...
    trans, r);

%-------------------------- END MAIN CODE -------------------------------%
end