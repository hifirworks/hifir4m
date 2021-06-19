function [x, flag, iter, ref_iter] = fgmresNull(A, b, ...
    us, leftnull, restart, rtol, maxit, M, ~, x0)
%fgmresNull Kernel for computing a null-space vector
%
% Syntax:
%   x = fgmresNull(A, b, us, leftnull, restart, rtol, maxit, M, ~, x0)
%   [x, flag, iter, ref_iter] = fgmresNull(...)
%
% See Also:
%   fgmres_HO

n = int32(size(b, 1));

% If RHS is zero, terminate
beta0 = sqrt(vec_sqnorm2(b));
if beta0 == 0
    x = zeros(n, 1);
    flag = int32(0);
    iter = int32(0);
    ref_iter = int32(0);
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

w = zeros(n, 1);

A_1nrm = norm(A,1);  % TODO: require a crs kernel for this

% condition number array
inrm_buf = zeros(restart+1, 1);

flag = int32(0);
iter = int32(0);
ref_iters = 0;
null_res_prev = realmax;
form_x_thres = min(sqrt(rtol),1e6);
for it_outer = 1:max_outer_iters
    % Compute the initial residual
    if it_outer > 1 || vec_sqnorm2(x) > 0
        w = ax_multiply(A, M, x, leftnull, w);
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
    % inner iteration steps
    N = bitshift(1, it_outer+3);
    inrm = 0;
    nrm = 0;
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
        %v = v/norm(v);

        % Store the preconditioned vector (A, M, N, u, b)
        [v, ref_iter, w] = iter_refine(A, M, N, v, w, leftnull, iter);
        ref_iters = ref_iters+ref_iter;
        % what we have G*q, z
        % NOTE, we also need to project of null space from preconditioned
        % subspace component
        v=mgs_null_proj(us,v);

        Z(:, j) = v;  % preconditioned subspace component
        w = ax_multiply(A, M, v, leftnull, w);

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
        % estimate the abs condition number of |R^{-T}|_\infty
        % TODO, change this to 2-norm
        [kappa, inrm, inrm_buf, nrm] = est_abs_cond(R, j, inrm, ...
            inrm_buf, nrm);

        if iter >= maxit
            flag = int32(1); % reached maxit
            break
        end
        iter = iter + 1;

        %m2c_printf('At iteration %d, |R^{-T}|_\\infty is %1.14e.\n', iter, resid);

        if j >= restart
            break;
        end

        if kappa >= form_x_thres
            x2 = x;
            y2 = backsolve(R, y, j);
            for i = 1:j
                x2 = x2 + y2(i) * Z(:, i);
            end
            x2=mgs_null_proj(us,x2);
            w = ax_multiply(A, M, x2, leftnull, w);
            er = vec_1norm(w)/(A_1nrm*vec_1norm(x2));
            %fprintf('|Ax|_1/|A|_1/|x|_1=%g\n', er);
            if er <= rtol
                flag=int32(2);
                break;
            end
            % safeguard preventing diverging
            if er >= null_res_prev
                j = j-1;
                flag = int32(2);
                ref_iters = ref_iters-ref_iter;
            end
            null_res_prev = er;
        end

        j = j + 1;
    end

    % Compute correction vector
    y = backsolve(R, y, j);
    for i = 1:j
        x = x + y(i) * Z(:, i);
    end
    
    % apply the projection for null space
    x=mgs_null_proj(us,x);

    if flag
        if flag == int32(2); flag = int32(0); end
        break;
    end
end

% Normalize the null space
x = x/sqrt(vec_sqnorm2(x));
ref_iter = int32(ref_iters);
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

function [kappa, inrm, inrm_buf, nrm] = est_abs_cond(R, i, ...
    inrm, inrm_buf, nrm)
coder.inline('always');
if i == 1
    inrm_buf(1) = 1./R(1,1);
    inrm = abs(inrm_buf(1));
    nrm = abs(R(1,1));
else
    if isempty(coder.target)
        s = inrm_buf(1:i-1)'*R(1:i-1,i);
    else
        s = 0.0;
        for j = int32(1):i-1
            s = s+inrm_buf(j)*R(j,i);
        end
    end
    k1 = 1-s; k2 = -1-s;
    if abs(k1)<abs(k2); inrm_buf(i) = k2/R(i,i); else; inrm_buf(i)=k1/R(i,i); end
    inrm = max(inrm, abs(inrm_buf(i)));
    nrm = max(nrm, sum(abs(R(1:i,i))));
end
kappa = inrm*nrm;
end

function [x, iter, w] = iter_refine(A, M, N, b, w, trans, gmres_iter)
coder.inline('always');
% applying iterative refinement A*x=b
% assuming M^{g} is a good approx of generalized inverse of A
if gmres_iter == 0
    N = min(N, 4);
end
bnorm = sqrt(vec_sqnorm2(b));
beta_L = 0.2;
beta_U = 10;
% beta_stag = 0;
% res_prev = 1;
x = zeros(size(b));
r = b;
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

function x = mgs_null_proj(us, x)
coder.inline('always');
if isempty(us); return; end
% Orthogonalize via MGS
for i = 1:size(us,2)
    u = us(:,i);
    x = x - (u'*x)*u;
end
end

function nrm1 = vec_1norm(v)
coder.inline('always');
n = int32(size(v, 1));
nrm1 = 0;
for i = int32(1):n; nrm1 = nrm1 + abs(v(i)); end
end
