function varargout = hilucsi4m_m_solve(dbase, varargin)
%HILUCSI4M_M_SOLVER - Accessing inv(M)
%
% Syntax:
%   x = hilucsi4m_m_solve(dbase, b)
%   x = hilucsi4m_m_solve(dbase, A, b, N)
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
%   x = hilucsi4m_m_solve(dbase, A, b, N) computes M\b in an stationary
%   iteration fashion, with inner steps of N.
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

if length(varargin) == 1
    [varargout{1:nargout}] = hilucsi4m_mex(HILUCSI4M_M_SOLVE, dbase, b);
else
    assert(length(varargin) == 3, 'inputs must be A, b, and N');
    b = varargin{2};
    N = varargin{3};
    if N < 2
        [varargout{1:nargout}] = hilucsi4m_mex(HILUCSI4M_M_SOLVE, dbase, b);
    else
        A = varargin{1};
        if issparse(A); A = hilucsi4m_sp2crs(A); end
        assert(isa(A.row_ptr, 'int32'));
        assert(isa(A.col_ind, 'int32'));
        [varargout{1:nargout}] = hilucsi4m_mex(HILUCSI4M_M_SOLVE, dbase, b, ...
            A.row_ptr, A.col_ind, A.val, N);
    end
end

%-------------------------- END MAIN CODE -------------------------------%
end