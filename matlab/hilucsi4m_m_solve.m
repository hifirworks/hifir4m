function varargout = hilucsi4m_m_solve(dbase, b)
%HILUCSI4M_M_SOLVER - Accessing inv(M)
%
% Syntax:
%   x = hilucsi4m_m_solve(dbase, b)
%   [x, t] = hilucsi4m_m_solve(___)
%
% Description:
%   HILUCSI4M_M_SOLVER allows one to directly access the preconditioner,
%   which is useful for developers who want to integrate HILUCSI into
%   existing (linear) solver frameworks.
%
%   x = hilucsi4m_m_solve(dbase, b) computes M\b given a factorized
%   database.
%
%   [x, t] = hilucsi4m_m_solve(___) can potentially report the
%   overhead-free wall-clock time.
%
% See Also:
%   HILUCSI4M_FACTORIZE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: GLPv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

[varargout{1:nargout}] = hilucsi4m_mex(HILUCSI4M_M_SOLVE, dbase, b);

%-------------------------- END MAIN CODE -------------------------------%
end