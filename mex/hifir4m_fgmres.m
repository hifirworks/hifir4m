function varargout = hifir4m_fgmres(dbase, A, b, varargin)
%HIFIR4M_FGMRES - right-preconditioned FGMRES with HIF
%
% Syntax:
%   x = hifir4m_fgmres(dbase, A, b)
%   x = hifir4m_fgmres(dbase, A, b, restart)
%   x = hifir4m_fgmres(dbase, A, b, restart, rtol)
%   x = hifir4m_fgmres(dbase, A, b, restart, rtol, maxit)
%   x = hifir4m_fgmres(dbase, A, b, restart, rtol, maxit, x0)
%   x = hifir4m_fgmres(dbase, A, b, restart, rtol, maxit, x0, verbose)
%   [x, flag] = hifir4m_fgmres(___)
%   [x, flag, iters] = hifir4m_fgmres(___)
%   [x, flag, iters, t] = hifir4m_fgmres(___)
%   [___] = hifir4m_fgmres(___, iter_refines)
%   [___] = hifir4m_fgmres(___, iter_refines, nsp_cst)
%
% Description:
%   HIFIR4M_GMRES is an optimized serial implementation of
%   right-preconditioned GMRES (sometimes referred as flexible GMRES [1]),
%   which has been shown to be more robust than the left-version [2]. The
%   GMRES here uses HIFIR as a right preconditioner, which works very
%   well (especially for saddle point and indifinite problems).
%
%   Before we dig into each of the syntaxes highlighted above, one thing
%   needs to be kept in mind is that the solver MUST be called only after
%   an instance of the preconditioner has been successfully factorized!
%
%   x = hifir4m_fgmres(dbase, A, b) computes the solution of A\b with
%   default parameters.
%
%   x = hifir4m_fgmres(dbase, A, b, restart, rtol, maxit) allows one to
%   customize restart (30), relative tolerance (1e-6) and maximum
%   iterations (500) for the GMRES solver; their default values are shown
%   in the parentheses.
%
%   x = hifir4m_fgmres(___, x0, verbose) further allows one to supply the
%   initial guess (all zeros) and verbose flag (true).
%
%   [x, flag, iters, t] = hifir4m_fgmres(___) indicates that there are
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
%   [___] = hifir4m_fgmres(___, iter_refines) indicates using iterative
%   refine kernel for the preconditioner.
%
%   [___] = hifir4m_fgmres(___, iter_refines, nsp_cst) solves a singular problem
%   with a (partial) constant mode that is specificed via a size-2 array
%   nsp_cst, in which the first entry is the starting const mode entry while
%   the ending index for the second element in nsp_cst
%
% Examples:
%   The following example shows how to use the GMRES solver
%       >> % assume we have dbase initialized and factorized
%       >> A = sprand(10, 10, 0.5);
%       >> b = rand(size(A, 1), 1);
%       >> x = hifir4m_fgmres(dbase, A, b);
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
%   HIFIR_FACTORIZE, GMRES

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: AGPLv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

gmres_pars = [30 1e-6 500]; % restart,rtol,maxit
x0 = [];
verbose = true;
iter_refines = int32(1);
for i = 1:min(3, length(varargin))
    if ~isempty(varargin{i}); gmres_pars(i) = varargin{i}; end
end
if length(varargin) > 3; x0 = varargin{4}; end
if length(varargin) > 4
    if ~isempty(varargin{5}); verbose = logical(varargin{5}); end
end
if length(varargin) > 5
    if ~isempty(varargin{6}); iter_refines = int32(varargin{6}); end
end
if isempty(x0)
    if isreal(b)
        x0 = zeros(size(b));
    else
        x0 = complex(zeros(size(b)), zeros(size(b)));
    end
end
% Convert to zero-based CRS
if issparse(A); A = hifir4m_sp2crs(A); end
A = hifir4m_ensure_int(A);
if length(varargin) < 7 || isempty(varargin{7})
    [varargout{1:nargout}] = hifir4m_mex(HIFIR4M_KSP_SOLVE, dbase, ...
        A.row_ptr, A.col_ind, A.val, b, gmres_pars(1), gmres_pars(2), ...
        gmres_pars(3), x0, verbose, iter_refines);
else
    [varargout{1:nargout}] = hifir4m_mex(HIFIR4M_KSP_SOLVE, dbase, ...
        A.row_ptr, A.col_ind, A.val, b, gmres_pars(1), gmres_pars(2), ...
        gmres_pars(3), x0, verbose, iter_refines, varargin{7});
end

%-------------------------- END MAIN CODE -------------------------------%
end