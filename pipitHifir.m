function [x, vs, us, flag, relres, iter, resids, times] = pipitHifir(A, b, ...
    nNull, varargin)
%pipitHifir - Compute PI solution of inconsistent systems.
%
% Syntax:
%   [x, vs, us] = pipitHifir(A, b, nNull, ...)
%   [x, vs, us, flag, relres, iter, resids, times] = pipitHifir(___)
%
% See Also:
%   pipitNull, fgmresNull

if nargin == 0
    help pipitHifir
    return;
end

p = inputParser;

% Initialize default arguments
addOptional(p, 'restart', int32(30), @(x) isempty(x) || isscalar(x) && (x>0));
addOptional(p, 'rtol', 1.e-6, @(x) isempty(x) || isscalar(x) && (x>0));
addOptional(p, 'maxit', int32(500), @(x) isempty(x) || isscalar(x) && (x>0));
addOptional(p, 'x0', cast([], class(b)), ...
    @(x) isempty(x) || isequal(size(x0), size(b)));
addOptional(p, 'rtol_null', 1e-12, @(x) isempty(x) || isscalar(x) && (x>0));
addParameter(p, 'verbose', int32(1), @isscalar);
addOptional(p, 'RNS', int32(1), @(x) isempty(x) || isscalar(x));

p.KeepUnmatched = true;
parse(p, varargin{:});

opts = p.Results;
% Process positional arguments
if isempty(opts.restart)
    opts.restart = int32(30);
end
if isempty(opts.maxit)
    opts.maxit = int32(500);
end
if isempty(opts.rtol)
    opts.rtol = 1.e-6;
end
if isempty(opts.rtol_null)
    opts.rtol_null = 1e-12;
end
if isempty(opts.RNS)
    opts.RNS = true;  % Whether or not skip computing RHS
end

% Create Hifir object
if opts.verbose
    fprintf(1, 'Computing hybrid factorization...\n');
end

times = zeros(4, 1);
args = namedargs2cell(p.Unmatched);
[hif, info, times(1)] = hifCreate(A, [], 'verbose', ...
    opts.verbose>1, args{:});

if opts.verbose
    fprintf(1, 'Finished ILU factorization in %.4g seconds \n', times(1));
    disp(info);
    fprintf(1, 'Starting Krylov solver for null space ...\n');
end

% Determine dimension of nullspace
if isempty(nNull) || nNull <= 0
    nNull = int32(1);
end

% Compute all left null space with a given known nullspace dimension
tic;
vs = pipitNull(A, hif, true, nNull, info.schur_rank ~= info.schur_size, ...
    opts.restart, opts.rtol_null, opts.maxit, opts.verbose);
times(2) = toc;

% Determine the least-square solution
b = mgs_null_proj(vs, b); % project off left nullspace

tic;
[x, flag, relres, iter, resids] = fgmres_MGS(A, b, opts.restart, ...
    opts.rtol, opts.maxit, hif, [], opts.x0);
times(3) = toc;

if opts.verbose
    if flag == 0
        fprintf(1, 'Finished solve LS in %d iterations and %.2f seconds.\n', iter, times(2));
    elseif flag == 3
        fprintf(1, 'GMRES stagnated after %d iterations and %.4g seconds.\n', iter, times(2));
    else
        fprintf(1, 'GMRES failed to converge after %d iterations and %.4g seconds.\n', iter, times(2));
    end
end

% Compute PI
us = [];
if ~opts.RNS; return; end
if issparse(A)
    % TODO: need to compute this for CRS
    if norm(A-A',1)/norm(A,1) <= 1e-12
        % symmetric
        us = vs;
    end
end

if isempty(us)
    % need to compute all the right null space
    tic;
    us = pipitNull(A, hif, false, size(vs, 2), ...
        info.schur_rank ~= info.schur_size, opts.restart, ...
        opts.rtol_null, opts.maxit, opts.verbose);
    times(4) = toc;
end

x = mgs_null_proj(us, x); % project off right nullspace

end

function x = mgs_null_proj(us, x)
% Orthogonalize via MGS
if isempty(us); return; end
for i = 1:size(us,2)
    u = us(:,i);
    x = x - (u'*x)*u;
end
end