function [vs, iters, ref_iters, times] = pipitNull(A, M, isleft, nNull, ...
    isSig, restart, rtol, maxit, verbose)
%pipitNull - Solving for nullspace
%
% Syntax:
%   vs = pipitNull(A, M, isleft, nNull, isSig, restart, rtol, maxit, verbose)
%   [vs, iters, ref_iters, times] = pipitNull(___)
%
% Description:
%   Small-dimension nullspace computation based on HIFIR and FGMRES. nNull
%   is the dimension of the nullspace. isSig indicates whether or not the
%   HIF preconditioner is singular, which can be determined by checking the
%   final Schur complement, i.e., info.schur_rank == info.schur_size.
%
% See Also:
%   fgmresNull

if issparse(A)
    n = int32(size(A, 1));
else
    n = A.nrows;
end
vs = create_rand_orth(n, nNull);
iters = zeros(nNull, 1, 'int32');
ref_iters = iters;
times = zeros(size(iters));
nullIter = int32(1);
while true
    if verbose
        fprintf(1, 'Looking for %d null-space vector ...\n', nullIter);
    end
    if ~isSig
        vs(:, nullIter) = iter_refine(A, M, 16, vs(:, nullIter), true);
        vs(:, nullIter) = vs(:, nullIter) / norm(vs(:, nullIter));
    end
    tn = tic;
    [vs(:, nullIter), flag, iters(nullIter), ref_iters(nullIter)] = ...
        fgmresNull(A, vs(:, nullIter), [], isleft, restart,  rtol, ...
        maxit, M, [], []);
    tn = toc(tn);
    times(nullIter) = tn;
    if verbose
        fprintf(1, 'Finished %d null-space vector (flag=%d) in %d(%d) iters and %.4gs.\n', ...
            nullIter, flag, iters(nullIter), ref_iters(nullIter), tn);
    end
    if flag
        nullIter = nullIter - 1;
        break;
    end
    if nullIter >= nNull; break; end
    nullIter = nullIter + 1;
end
if nullIter <= 0; vs = []; end
if ~isempty(vs)
    if nullIter < size(vs, 2)
        vs = vs(:, 1:nullIter);
    end
    [vs, ~] = qr(vs, 0);  % reorth
end
end

function vs = create_rand_orth(n, m)
[vs, ~] = qr(rand(n, m), 0);
end

function b = ax_multiply(A, M, x, trans, b)
coder.inline('always');
if issparse(A)
    if ~trans; b = A * x; else; b = A' * x; end
else
    if ~trans
        b = crs_prodAx(M.A, x, b);
    else
        b = crs_prodAtx(M.A, x, b);
    end
end
end

function [x, iter] = iter_refine(A, M, N, b, trans)
% applying iterative refinement A*x=b
bnorm = sqrt(vec_sqnorm2(b));
beta_L = 0.2;
beta_U = 1e8;
% beta_stag = 0;
% res_prev = 1;
x = zeros(size(b));
r = b;
w = b;
for iter = int32(1):N
    if isa(M, 'function_handle')
        % TODO: handle transpose
        r = M(r);
    else
        if trans; op = 'SH'; else; op = 'S'; end
        r = M.apply(r, op, -1);
    end
    x = x+r;
    w = ax_multiply(A, M, x, trans, w);
    r = b - w;
    res = sqrt(vec_sqnorm2(r))/bnorm;
    if res > beta_U
        break;
    end
    if res <= beta_L; break; end
end
end