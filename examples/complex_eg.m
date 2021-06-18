%{
    This example code uses HIF-preconditioned GMRES(30) to solve a
    complex testing system "kim1" from SuiteSparse Matrix Collection.
%}

clear;

load('kim1.mat', 'Problem');

% RHS is A*1
A = Problem.A;
n = size(A, 1);
b = A*ones(n, 1);

% solve with HIF-preconditioned GMRES
[x, flag] = gmresHif(A, b, 'is_complex', true);
if flag
    fprintf(2, 'warning! solver failed with flag %d\n', flag);
end
res = norm(b - A * x) / norm(b);
if res > 1e-6
    fprintf(2, 'residual too large %.4g\n', res);
end
