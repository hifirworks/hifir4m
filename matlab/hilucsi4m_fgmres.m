function varargout = hilucsi4m_fgmres(dbase, A, b, varargin)
%HILUCSI4M_FGMRES - Flexible (right-preconditioned) GMRES with HILUCSI
%
% Syntax:
%   x = hilucsi4m_fgmres(dbase, A, b)
%   x = hilucsi4m_fgmres(dbase, A, b, restart)
%   x = hilucsi4m_fgmres(dbase, A, b, restart, rtol)
%   x = hilucsi4m_fgmres(dbase, A, b, restart, rtol, maxit)
%   x = hilucsi4m_fgmres(dbase, A, b, restart, rtol, maxit, x0)
%   x = hilucsi4m_fgmres(dbase, A, b, restart, rtol, maxit, x0, verbose)
%   [x, flag] = hilucsi4m_fgmres(___)
%   [x, flag, iters] = hilucsi4m_fgmres(___)
%   [x, flag, iters, t] = hilucsi4m_fgmres(___)
%   [___] = hilucsi4m_fgmres(___, update)
%   [___] = hilucsi4m_fgmres(___, update, nsp_cst)
%
% Description:
%   HILUCSI4M_FGMRES is an optimized serial implementation of
%   right-preconditioned GMRES (sometimes referred as flexible GMRES [1]),
%   which has been shown to be more robust than the left-version [2]. The
%   FGMRES here uses HILUCSI as a right preconditioner, which works very
%   well (especially for saddle point and indifinite problems).
%
%   Before we dig into each of the syntaxes highlighted above, one thing
%   needs to be kept in mind is that the solver MUST be called only after
%   an instance of the preconditioner has been successfully factorized!
%
%   x = hilucsi4m_fgmres(dbase, A, b) computes the solution of A\b with
%   default parameters.
%
%   x = hilucsi4m_fgmres(dbase, A, b, restart, rtol, maxit) allows one to
%   customize restart (30), relative tolerance (1e-6) and maximum
%   iterations (500) for the GMRES solver; their default values are shown
%   in the parentheses.
%
%   x = hilucsi4m_fgmres(___, x0, verbose) further allows one to supply the
%   initial guess (all zeros) and verbose flag (true).
%
%   [x, flag, iters, t] = hilucsi4m_fgmres(___) indicates that there are
%   three more optoinal outputs, namely
%       flag - solver status flag
%           0  - successed
%           1  - diverged
%           2  - stagnated
%           3  - break-down
%           <0 - input errors
%       iters - total iterations used
%       t - solving wall-clock time (no MATLAB interperater overhead)
%
%   [___] = hilucsi4m_fgmres(___, update) indicates using updated kernel
%   for the preconditioner.
%
%   [___] = hilucsi4m_fgmres(___, update, nsp_cst) solves a singular problem
%   with a (partial) constant mode that is specificed via a size-2 array
%   nsp_cst, in which the first entry is the starting const mode entry while
%   the ending index for the second element in nsp_cst
%
% Examples:
%   The following example shows how to use the FGMRES solver
%       >> % assume we have dbase initialized and factorized
%       >> A = sprand(10, 10, 0.5);
%       >> b = rand(size(A, 1), 1);
%       >> x = hilucsi4m_fgmres(dbase, A, b);
%       >> assert(norm(A-b)/norm(b)<=1e-6);
%
% References:
%   [1] Saad, Y. (1993). A flexible inner-outer preconditioned GMRES
%       algorithm. SIAM Journal on Scientific Computing, 14(2), 461-469.
%   [2] Ghai, A., Lu, C., & Jiao, X. (2019). A comparison of preconditioned
%       Krylov subspace methods for large‐scale nonsymmetric linear systems.
%       Numerical Linear Algebra with Applications, 26(1), e2215.
%
% See Also:
%   HILUCSI_FACTORIZE, GMRES

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: GLPv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

gmres_pars = [30 1e-6 500]; % restart,rtol,maxit
x0 = [];
verbose = true;
update = false;
for i = 1:min(3, length(varargin))
    if ~isempty(varargin{i}); gmres_pars(i) = varargin{i}; end
end
if length(varargin) > 3; x0 = varargin{4}; end
if length(varargin) > 4
    if ~isempty(varargin{5}); verbose = logical(varargin{5}); end
end
if length(varargin) > 5
    if ~isempty(varargin{6}); update = logical(varargin{6}); end
end
if isempty(x0); x0 = zeros(size(b)); end
% Convert to zero-based CRS
if issparse(A); A = hilucsi4m_sp2crs(A); end
assert(isa(A.row_ptr, 'int32'));
assert(isa(A.col_ind, 'int32'));
A = hilucsi4m_2int64(A);
if length(varargin) < 7 || isempty(varargin{7})
    [varargout{1:nargout}] = hilucsi4m_mex(HILUCSI4M_KSP_SOLVE, dbase, ...
        A.row_ptr, A.col_ind, A.val, b, gmres_pars(1), gmres_pars(2), ...
        gmres_pars(3), x0, verbose, update);
else
    [varargout{1:nargout}] = hilucsi4m_mex(HILUCSI4M_KSP_SOLVE, dbase, ...
        A.row_ptr, A.col_ind, A.val, b, gmres_pars(1), gmres_pars(2), ...
        gmres_pars(3), x0, verbose, update, varargin{7});
end

%-------------------------- END MAIN CODE -------------------------------%
end