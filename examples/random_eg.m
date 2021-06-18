%{
    This example uses a density 0.5 of 10x10 sparse matrix to show how
    to use HIF-preconditioned GMRES solver.
%}
clear;

% setup testing matrix
A = sprand(10, 10, 0.5);
b = rand(10, 1);

% Call HIF-preconditioned GMRES
[x, flag] = gmresHif(A, b, 'is_mixed', true);  % using mixed-precision
if flag
    fprintf(2, 'warning! solver failed with flag %d\n', flag);
end
res = norm(b - A * x) / norm(b);
if res > 1e-6
    fprintf(2, 'residual too large %.4g\n', res);
end
    