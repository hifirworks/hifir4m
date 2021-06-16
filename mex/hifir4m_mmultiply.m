function varargout = hifir4m_mmultiply(dbase, varargin)
%HIFIR4M_MMULTIPLY - Accessing inv(M)
%
% Syntax:
%   x = hifir4m_mmultiply(dbase, b)
%   x = hifir4m_mmultiply(dbase, b, trans)
%   x = hifir4m_mmultiply(dbase, b, trans, r)
%   [x, t] = hifir4m_mmultiply(___)
%
% Description:
%   HIFIR4M_MMULTIPLY computes the multilevel matrix-vector multiplication.
%
%   x = hifir4m_mmultiply(dbase, b) computes M*b given a factorized
%   database.
%
%   x = hifir4m_mmultiply(dbase, b, trans) computes M'*b given a factorized
%   database if trans==true.
%
%   x = hifir4m_mmultiply(dbase, b, trans, r) computes a multilevel matrix-
%   vector multiplication with operation trans with a customized final Schur
%   complement dimension r.
%
%   [x, t] = hifir4m_mmultiply(___) can potentially report the
%   overhead-free wall-clock time.
%
% See Also:
%   HIFIR4M_SOLVE

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
[varargout{1:nargout}] = hifir4m_mex(HIFIR4M_M_MULTIPLY, dbase, varargin{1}, ...
    trans, r);

%-------------------------- END MAIN CODE -------------------------------%
end