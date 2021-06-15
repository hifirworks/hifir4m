function [x, flag, relres, iter, resids] = fgmres_HO(A, b, ...
    restart, rtol, maxit, M, ~, x0)
%fgmres_HO   Kernel of FGMRES using Householder algorithm
%
%   x = fgmres_HO(A, b, restart, rtol, maxit, M, ~, x0)
%
%   [x, flag, iter, resids] = fgmres_HO(...)
%
% See also: gmres, fgmres_CGS, fgmres_MGS

% Note: The algorithm uses Householder reflectors for orthogonalization.
% It is more expensive than Gram-Schmidt but is more robust.

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

% Householder matrix or upper-triangular matix
V = zeros(n, restart);
R = zeros(restart, restart);

% Temporary solution
y = zeros(restart+1, 1);

% Preconditioned subspace
Z = zeros(n, restart);

% Given's rotation vectors
J = zeros(2, restart);

% Corrections at outer loop
if nargout > 4
    resids = zeros(maxit, 1);
end

w = zeros(n, 1);

flag = int32(0);
iter = int32(0);
relres = 1;
for it_outer = 1:max_outer_iters
    % Compute the initial residual
    if it_outer > 1 || vec_sqnorm2(x) > 0
        if issparse(A)
            w = A * x;
        else
            w = crs_prodAx(M.A, x, w);
        end
        u = b - w;
    else
        u = b;
    end

    beta2 = vec_sqnorm2(u);
    beta = sqrt(beta2);

    % Prepare the first Householder vector
    if u(1) < 0
        beta = -beta;
    end
    updated_norm = sqrt(2*beta2+2*real(u(1))*beta);
    u(1) = u(1) + beta;
    u = u / updated_norm;

    % The first Householder entry
    y(1) = - beta;
    V(:, 1) = u;

    j = int32(1);
    while true
        % Construct the last vector from the Householder reflectors

        %  v = Pj*ej = ej - 2*u*u'*ej
        v = -2 * conj(V(j, j)) * V(:, j);
        v(j) = v(j) + 1;
        %  v = P1*P2*...Pjm1*(Pj*ej)
        if isempty(coder.target)
            % This is faster when interpreted
            for i = (j - 1): - 1:1
                v = v - 2 * (V(:, i)' * v) * V(:, i);
            end
        else
            % This is faster when compiled
            for i = (j - 1): - 1:1
                s = conj(V(i, i)) * v(i);
                for k = i + 1:n
                    s = s + conj(V(k, i)) * v(k);
                end
                s = 2 * s;

                for k = i:n
                    v(k) = v(k) - s * V(k, i);
                end
            end
        end
        %  Explicitly normalize v to reduce the effects of round-off.
        v = v / sqrt(vec_sqnorm2(v));

        % Compute the preconditioned vector
        if isequal(class(M), 'function_handle')
            v = M(v);
        else
            v = M.apply(v);
        end

        Z(:, j) = v;
        if issparse(A)
            w = A * v;
        else
            w = crs_prodAx(M.A, v, w);
        end

        % Orthogonalize the Krylov vector
        %  Form Pj*Pj-1*...P1*Av.
        if isempty(coder.target)
            % This is faster when interpreted
            for i = 1:j
                w = w - 2 * (V(:, i)' * w) * V(:, i);
            end
        else
            % This is faster when compiled
            for i = 1:j
                s = conj(V(i, i)) * w(i);
                for k = i + 1:n
                    s = s + conj(V(k, i)) * w(k);
                end
                s = s * 2;

                for k = i:n
                    w(k) = w(k) - s * V(k, i);
                end
            end
        end

        % Update the rotators
        % Determine Pj+1.
        if j < n
            %  Construct u for Householder reflector Pj+1.
            u(j) = 0;
            u(j+1) = w(j+1);
            alpha2 = conj(w(j+1)) * w(j+1);
            for k = j + 2:n
                u(k) = w(k);
                alpha2 = alpha2 + conj(w(k)) * w(k);
            end

            if alpha2 > 0
                alpha = sqrt(alpha2);
                if u(j+1) < 0
                    alpha = -alpha;
                end
                if j < restart
                    updated_norm = sqrt(2*alpha2+2*real(u(j+1))*alpha);
                    u(j+1) = u(j+1) + alpha;
                    for k = j + 1:n
                        V(k, j+1) = u(k) / updated_norm;
                    end
                end

                %  Apply Pj+1 to v.
                w(j+2:end) = 0;
                w(j+1) = - alpha;
            end
        end

        %  Apply Given's rotations to the newly formed v.
        for colJ = 1:j - 1
            tmpv = w(colJ);
            w(colJ) = conj(J(1, colJ)) * w(colJ) + conj(J(2, colJ)) * w(colJ+1);
            w(colJ+1) = - J(2, colJ) * tmpv + J(1, colJ) * w(colJ+1);
        end

        %  Compute Given's rotation Jm.
        if j < n
            rho = sqrt(w(j)'*w(j)+w(j+1)'*w(j+1));
            J(1, j) = w(j) ./ rho;
            J(2, j) = w(j+1) ./ rho;
            y(j+1) = - J(2, j) .* y(j);
            y(j) = conj(J(1, j)) .* y(j);
            w(j) = rho;
        end

        R(1:j, j) = w(1:j);

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
