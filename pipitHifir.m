function [x, us, vs, flag, relres, iter, resids, times] = ...
    pipitHifir(A, b, nNull, varargin)
%pipitHifir - Compute pseudoinverse solution of an inconsistent system.
%
%      x = pipitHifir(A, b, nNull, [, restart, rtols, maxit, x0, vs])
%    solves an inconsistent sparse linear system using PIPIT with
%    HIFIR as right preconditioner. Matrix A can be in MATLAB's built-in
%    sparse format or in CRS struct with row_ptr, col_ind, and vals fields.
%
%    nNull specifies the dimension of the null-space. (This argument is
%        presently required but will be optional in the future.)
%
%    The optional arguments include:
%    M: user-provided HIF preconditioner. The default value is empty indicating
%        that this function will compute a preconditioner.
%
%    restart: the number of iterations before GMRES restarts. If empty or
%        missing, the default value is 30.
%
%    rtols: a two-vector specifying the relative tolerances for the
%        solution and nulll-space vectors. If empty or missing, the default
%        values are [1.e-6, 1.e-12].
%
%    maxit: the maximum number of iterations. If empty or missing,
%       the default value is 500. (Note that unlike MATLAB's built-in gmres,
%       maxit here refers to the total number of iterations.)
%
%    x0: the initial guess.
%    vs: the right null space, if known. If vs is specified but is a scalar
%       zero, PIPIT will compute the least-squares solution (i.e., it will
%       skip projecting off the right-null-space component in the solution).
%
%       x = pipitHifir(A, b, nNull, ..., 'name', value, ...)
%    allows omitting some of the optional arguments followed by name-value
%    pairs for the parameters. The parameters include:
%       verbose:  verbose level. Default is 1.
%    Unmatched parameters are passed to HifParams.
% 
%    [x, us, vs, flag, relres, iter, resids, times] = pipitHifir(A, b, nNull, ...)
%    returns additional output in addition to the solution vector.
%
%    us: the (computed) left null-space basis vectors
%    vs: the (computed) right null-space basis vectors
%    flag: 0 - converged to the desired tolerance TOL within MAXIT iterations.
%          1 - iterated maxit times but did not converge.
%          3 - stagnated (two consecutive iterates were the same).
%    relres: the relative residual.
%    iter:   the number of iterations.
%    resids: contains the residial history of all iterations.
%    times:  returns the setup time and solve time in seconds in a 4-vector.
%
% See also gmresHif, fgmresNull, orthNulls

if nargin == 0
    help pipitHifir
    return;
end

p = inputParser;

% Initialize default arguments
addRequired(p, 'nNull', @(x) isscalar(x) && (x>0));
addParameter(p, 'M', [], @(x) isempty(x) || isa(x, 'Hifir'));
addParameter(p, 'restart', int32(30), @(x) isempty(x) || isscalar(x));
addParameter(p, 'rtols', [1.e-6,1.e-12], @(x) numel(x)<=2);
addParameter(p, 'maxit', int32(500), @(x) isempty(x) || isscalar(x));
addParameter(p, 'x0', cast([], class(b)), ...
    @(x) isempty(x) || isequal(size(x), size(b)));
addParameter(p, 'vs', cast([], class(b)), ...
    @(x) isempty(x) || isscalar(x) && x==0 || isequal(size(x), size(b)));
addParameter(p, 'verbose', int32(1), @isscalar);

p.KeepUnmatched = true;
parse(p, nNull, varargin{:});

opts = p.Results;
% Process positional arguments
if isempty(opts.restart) || opts.restart <= 0
    opts.restart = int32(30);
end
if isempty(opts.maxit) || opts.maxit <= 0
    opts.maxit = int32(500);
end
if isempty(opts.rtols) || opts.rtols(1) <= 0
    rtol = 1.e-6;
else
    rtol = opts.rtols(1);
end
if isempty(opts.rtols) || numel(opts.rtols) < 2 || opts.rtols(2) <= 0
    rtol_null = 1.e-12;
else
    rtol_null = opts.rtols(2);
end

% Create Hifir object
computedHif = isempty(opts.M);
times = zeros(4, 1);

if computedHif
    if opts.verbose
        fprintf(1, 'Computing hybrid incomplete factorization...\n');
    end
    args = namedargs2cell(p.Unmatched);
    [hif, info, times(1)] = hifCreate(A, [], 'verbose', opts.verbose>1, args{:});
else
    if opts.verbose
        fprintf(1, 'Using user-provided hybrid incomplete factorization...\n');
    end
    hif = opts.M;
    % get info
    info = hifir4m_mex(HifEnum.QUERY, hif.hdl);
end
if opts.verbose
    if computedHif
        fprintf(1, 'Finished in %.4g seconds \n', times(1));
        disp(info);
    end
    fprintf(1, 'Computing left null space ...\n');
end

% Determine the dimension of nullspace
if isempty(nNull) || nNull <= 0
    nNull = int32(1);
end

% Compute the left null space using the given nullspace dimension
tic;
[us, iters_left, iters_ir_left] = orthNulls(A, hif, nNull, true, ...
    info.schur_rank ~= info.schur_size, opts.restart, rtol_null, opts.maxit);
b = mgs_null_proj(us, b); % project off left nullspace
times(2) = toc;

% Compute the least-square solution
if opts.verbose
    fprintf(1, 'Finished after %d FGMRES iterations with %d iterative refinements in %.4g seconds.\n', ...
        sum(iters_left), sum(iters_ir_left), times(2));
    fprintf(1, 'Computing least-squares solution...\n');
end

tic;
[x, flag, relres, iter, resids] = fgmres_MGS(A, b, opts.restart, ...
    rtol, opts.maxit, hif, [], opts.x0);
times(3) = toc;

if opts.verbose
    if flag == 0
        fprintf(1, 'Finished after %d GMRES iterations in %.4g seconds.\n', iter, times(3));
    elseif flag == 3
        fprintf(1, 'GMRES stagnated after %d iterations and %.4g seconds.\n', iter, times(3));
    else
        fprintf(1, 'GMRES failed to converge after %d iterations and %.4g seconds.\n', iter, times(3));
    end
end

% Return least-squares solution if opts.vs is 0
if isscalar(opts.vs)
    vs = [];
    return;
end

% Compute pseudo-inverse solution
iters_right = [];
if isempty(opts.vs)
    if issparse(A)
        if norm(A-A',1) <= rtol_null * norm(A,1)
            if opts.verbose
                if isreal(A)
                    fprintf(1, 'Matrix is symmetric.');
                else
                    fprintf(1, 'Matrix is Hermitian.');
                end
                fprintf(1, ' Reusing left null-space for right null space.\n');
            end
            tic;
            vs = us;
        else
            % compute the right null space
            if opts.verbose
                fprintf(1, 'Converting least-squares solution to pseudo-inverse solution...\n');
            end
            tic;
            [vs, iters_right, iters_ir_right] = orthNulls(A, hif, size(us, 2), false, ...
                info.schur_rank ~= info.schur_size, opts.restart, ...
                rtol_null, opts.maxit);
        end
    else
        % TODO: need to implement this for CRS
        assert(false, 'CRS is not yet implemented');
    end
else
    vs = opts.vs;
end

x = mgs_null_proj(vs, x); % project off right nullspace
times(4) = toc;
if opts.verbose
    if ~isempty(iters_right)
        fprintf(1, 'Finished after %d FGMRES iterations with %d iterative refinements in %.4g seconds.\n', ...
            sum(iters_right), sum(iters_ir_right), times(4));
    else
        fprintf(1, 'Finished in %g seconds.\n', times(4));
    end
end

end

function x = mgs_null_proj(vs, x)
% Orthogonalize via MGS
if isempty(vs); return; end
for i = 1:size(vs,2)
    u = vs(:,i);
    x = x - (u'*x)*u;
end
end

function test %#ok<DEFNU>
%!test
%! A = gallery('neumann', 256);
%! b = rand(size(A, 1), 1);
%! x = pipitHifir(A, b, 1);
%! assert(norm(A'*(b-A*x))/norm(A'*b) <= 1e8);

end
