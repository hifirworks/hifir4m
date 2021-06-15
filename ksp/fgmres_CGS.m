function [x, flag, relres, iter, resids] = fgmres_CGS(A, b, ...
    restart, rtol, maxit, M, ~, x0)
% fgmres_CGS    Kernel of FGMRES using classical Gram-Schmidt
%
%   x = fgmres_CGS(A, b, restart, rtol, maxit, M, ~, x0)
%
%   [x, flag, iter, resids] = fgmres_CGS(...)
%
% See also: gmres, fgmres_MGS, fgmres_HO

% Note: The algorithm uses the classical Gram-Schmidt orthogonalization.
% It has more parallelism than modified  Gram-Schmidt but is less stable.
% It is also less stable than the Householder algorithm.

n = int32(size(b, 1));

% If RHS is zero, terminate
beta0 = sqrt(vec_sqnorm2(b));
if beta0 == 0
    x = zeros(n, 1);
    flag = int32(0);
    iter = int32(0);
    relres = 0;
    resids = 0;
    return;
end

% Number of inner iterations
if restart > n
    restart = n;
elseif restart <= 0
    restart = int32(1);
end

% Determine the maximum number of outer iterations
max_outer_iters = int32(ceil(double(maxit)/double(restart)));

% Initialize x
if isempty(x0)
    x = zeros(n, 1);
else
    x = x0;
end

% Local linear system
y = zeros(restart+1, 1);
R = zeros(restart, restart);

% Orthognalized Krylov subspace
Q = zeros(n, restart);

% Preconditioned subspace
Z = zeros(n, restart);

% Given's rotation vectors
J = zeros(2, restart);

% Buffer spaces
v = zeros(n, 1);

if nargout > 4
    resids = zeros(maxit, 1);
end

flag = int32(0);
iter = int32(0);
relres = 1;
for it_outer = 1:max_outer_iters
    % Compute the initial residual
    if it_outer > 1 || vec_sqnorm2(x) > 0
        if issparse(A)
            v = A * x;
        else
            v = crs_prodAx(M.A, x, v);
        end
        v = b - v;
    else
        v = b;
    end

    beta2 = vec_sqnorm2(v);
    beta = sqrt(beta2);

    % The first Q vector
    y(1) = beta;
    Q(:, 1) = v / beta;

    j = int32(1);
    while true
        w = Q(:, j);
        % Compute the preconditioned vector
        if isequal(class(M), 'function_handle')
            w = M(w);
        else
            w = M.apply(w);
        end

        % Store the preconditioned vector
        Z(:, j) = w;
        if issparse(A)
            v = A * w;
        else
            v = crs_prodAx(M.A, w, v);
        end

        % Perform classical Gram-Schmidt orthogonalization
        w = v;
        for k = 1:j
            R(k, j) = w' * Q(:, k);
            v = v - R(k, j) * Q(:, k);
        end

        vnorm2 = vec_sqnorm2(v);
        vnorm = sqrt(vnorm2);
        if j < restart
            Q(:, j+1) = v / vnorm;
        end

        %  Apply Given's rotations to R(:,j)
        for colJ = 1:j-1
            tmpv = R(colJ, j);
            R(colJ, j) = conj(J(1, colJ)) * R(colJ, j) + conj(J(2, colJ)) * R(colJ+1, j);
            R(colJ+1, j) = - J(2, colJ) * tmpv + J(1, colJ) * R(colJ+1, j);
        end

        %  Compute Given's rotation Jm.
        rho = sqrt(R(j, j)'*R(j, j)+vnorm2);
        J(1, j) = R(j, j) ./ rho;
        J(2, j) = vnorm ./ rho;
        y(j+1) = - J(2, j) .* y(j);
        y(j) = conj(J(1, j)) .* y(j);
        R(j, j) = rho;

        resid_prev = relres;
        relres = abs(y(j+1)) / beta0;
        if relres >= resid_prev * (1 - 1.e-8)
            flag = int32(3); % stagnated
            break
        elseif iter >= maxit
            flag = int32(1); % reached maxit
            break
        end
        iter = iter + 1;

        % save the residual
        if nargout > 4
            resids(iter) = relres;
        end

        if relres < rtol || j >= restart
            break;
        end
        j = j + 1;
    end

    % Compute correction vector
    y = backsolve(R, y, j);
    for i = 1:j
        x = x + y(i) * Z(:, i);
    end

    if relres < rtol || flag
        break;
    end
end

if nargout > 4
    resids = resids(1:iter);
end

if relres <= rtol * (1 + 1.e-8)
    flag = int32(0);
end

end
