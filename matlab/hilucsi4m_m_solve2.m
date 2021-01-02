function [x, t] = hilucsi4m_m_solve2(dbase, b)
%HILUCSI4M_M_SOLVER2 - Accessing inv(M) with two RHS
%
% Syntax:
%   x = hilucsi4m_m_solve2(dbase, b)
%   [x, t] = hilucsi4m_m_solve(___)
%
% Description:
%   HILUCSI4M_M_SOLVER2 allows one to directly access the preconditioner,
%   which is useful for developers who want to integrate HILUCSI into
%   existing (linear) solver frameworks. Note that this routines works with
%   two RHS.
%
%   x = hilucsi4m_m_solve(dbase, b) computes M\b given a factorized
%   database.
%
%   [x, t] = hilucsi4m_m_solve(___) can potentially report the
%   overhead-free wall-clock time.
%
% See Also:
%   HILUCSI4M_M_SOLVE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: GLPv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

assert(size(b, 2) == 2, 'RHS must be n-by-2');
% transpose for better cache performance
b = b.';

[x, t] = hilucsi4m_mex(HILUCSI4M_M_SOLVE2, dbase, b);
x = x.';

%-------------------------- END MAIN CODE -------------------------------%
end