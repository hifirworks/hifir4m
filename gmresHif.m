function [x, flag, relres, iter, resids, times] = gmresHif(A, b, varargin)
% gmresHif  Restarted GMRES with HIF as right preconditioner
%
%       x = gmresHif(A, b [, M, restart, rtol, maxit, x0])
%
%    solves a sparse linear system using HIF as the preconditioner.
%    Matrix A can be in MATLAB's built-in sparse format or in CRS struct
%    with row_ptr, col_ind, and vals fields.
%
%    restart specifies the number of iterations before GMRES restarts.
%    If empty or missing, the default value is 30.
%
%    rtol specifies the relative tolerance. If empty or missing,
%    the default value is 1.e-6.
%
%    maxit specifies the maximum number of iterations. If empty or missing,
%    the default value is 500. (Note that unlike MATLAB's built-in gmres,
%    maxit here refers to the total number of iterations.)
%
%    x0 specifies the initial guess.
%
%       x = gmresHif(A, b, ..., 'name', value, ...)
%    allows omitting some of the optional arguments followed by name-value
%    pairs for the parameters. The parameters include:
%       M:        HIF preconditioner. Default is [], which means this
%                 function will compute a HIF preconditioner. For dynamic
%                 and/or nonlinear problem, it is preferable to allow
%                 the users to manage preconditioners outside of this
%                 function.
%       verbose:  verbose level. Default is 1.
%       pcside:   Side of preconditioner ('left' or 'right'). If 'left',
%                 MATLAB built-in GMRES will be called. Default is 'right'.
%       orth:     Orthogonalization method ('MGS', 'CGS', or 'HO') for
%                 right-preconditioned GMRES. Default is 'MGS'.
%    Unmatched parameters are passed to HifParams.
% 
%       [x, flag, relres, iter, resds, times] = gmresHif(...)
%    returns a convergence flag.
%
%    flag: 0 - converged to the desired tolerance TOL within MAXIT iterations.
%          1 - iterated maxit times but did not converge.
%          3 - stagnated (two consecutive iterates were the same).
%    relres: the relative residual.
%    iter:   the number of iterations.
%    resids: contains the residial history of all iterations.
%    times:  returns the setup time and solve time in seconds in a 2-vector.

if nargin == 0
    help gmresHif
    return;
end

p = inputParser;

% Initialize default arguments
addParameter(p, 'M', [], @(x) isempty(x) || isa(x, 'Hifir'));
addParameter(p, 'restart', int32(30), @(x) isempty(x) || isscalar(x));
addParameter(p, 'rtol', 1.e-6, @(x) isempty(x) || isscalar(x));
addParameter(p, 'maxit', int32(500), @(x) isempty(x) || isscalar(x));
addParameter(p, 'x0', cast([], class(b)), ...
    @(x) isempty(x) || isequal(size(x), size(b)));

addParameter(p, 'verbose', int32(1), @isscalar);
addParameter(p, 'pcside', 'right', ...
    @(position) strcmp(position, 'left') || strcmp(position, 'right'));
addParameter(p, 'orth', 'MGS', ...
    @(x) isequal(x, 'MGS') || isequal(x, 'CGS') || isequal(x, 'HO'));

p.KeepUnmatched = true;
parse(p, varargin{:});

opts = p.Results;
% Process positional arguments
if isempty(opts.restart) || opts.restart <= 0
    opts.restart = 30;
end
if isempty(opts.maxit) || opts.maxit <= 0
    opts.maxit = 500;
end
if isempty(opts.rtol) || opts.rtol <= 0
    opts.rtol = 1.e-6;
end

% Create Hifir object
computedHif = isempty(opts.M);
times = zeros(2, 1);
if computedHif
    if opts.verbose
        fprintf(1, 'Computing hybrid incomplete factorization...\n');
    end
    args = namedargs2cell(p.Unmatched);
    [hif, info, times(1)] = hifCreate(A, [], 'verbose', ...
        opts.verbose>1, args{:});
else
    if opts.verbose
        fprintf(1, 'Using user-provided hybrid incomplete factorization...\n');
    end
    hif = opts.M;
end
M = @(x) hifApply(hif, x);

if opts.verbose
    if computedHif
        fprintf(1, 'Computed HIF factorization in %.4g seconds \n', times(1));
        disp(info);
    end
    fprintf(1, 'Starting GMRES with HIF as right-preconditioner...\n');
end

tic;
if strcmp(opts.pcside, 'left')
    max_outer_iters = int32(ceil(double(opts.maxit)/double(opts.restart)));
    [x, flag, relres, iters, resids] = gmres(A, b, ...
        opts.restart, opts.rtol, max_outer_iters, M, [], opts.x0);
    iter = (iters(1)-1)*opts.restart+iters(2);
else
    kernel = ['fgmres_', opts.orth];
    kernel_func = eval(['@' kernel]);

    [x, flag, relres, iter, resids] = kernel_func(A, b, ...
        opts.restart, opts.rtol, opts.maxit, M, [], opts.x0);
end
times(2) = toc;

if opts.verbose
    if flag == 0
        fprintf(1, 'Computed solution in %d iterations and %.4g seconds.\n', iter, times(2));
    elseif flag == 3
        fprintf(1, 'GMRES stagnated after %d iterations and %.4g seconds.\n', iter, times(2));
    else
        fprintf(1, 'GMRES failed to converge after %d iterations and %.4g seconds.\n', iter, times(2));
    end
end

end

function test
%!test
%!shared A, b, rtol
%! % system('gd-get -O -p 0ByTwsK5_Tl_PemN0QVlYem11Y00 fem2d"*".mat');
%! s = load('fem2d_cd.mat');
%! A = s.A;
%! s = load('fem2d_vec_cd.mat');
%! b = s.b;
%! rtol = 1.e-5;
%
%! [x, flag, relres, iter, resids] = gmresHif(A, b, [], rtol, 100, 'pcside', 'left');
%! assert(norm(b - A*x) <= rtol * norm(b))
%!
%!test
%! [x, flag, relres, iter, resids] = gmresHif(A, b, [], rtol, 100, 'orth', 'MGS');
%! assert(norm(b - A*x) <= rtol * norm(b))

%!test
%! [x, flag, relres, iter, resids] = gmresHif(A, b, [], rtol, 100, 'orth', 'CGS');
%! assert(norm(b - A*x) <= rtol * norm(b))

%!test
%! [x, flag, relres, iter, resids] = gmresHif(A, b, [], rtol, 100, 'orth', 'HO');
%! assert(norm(b - A*x) <= rtol * norm(b))

end
