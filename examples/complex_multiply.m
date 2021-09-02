%{
    This example code demonstrates how to apply multilevel triangular
    solve and matrix-vector multiplication in both regular and Hermitian
    modes for the testing matrix "kim1" from SuiteSparse Matrix
    Collections.
%}

clear;

load('kim1.mat', 'Problem');

A = Problem.A;
n = size(A,1);
b = A*ones(n,1);

% Compute HIF preconditioner
[hdl, stats] = hifCreate(A, [], 'is_complex', true);
fprintf(1, 'Factorization finished with stats:\n\n');
disp(stats);

% Apply solve and multiplication
x = hifApply(hdl, b, 'S');  % triangular solve
b2 = hifApply(hdl, x, 'M'); % mat-vec multiplication
res = norm(b2 - b) / norm(b);
assert(res <= 1e-10, 'residual %g is too large', res);

% Apply solve and multiplication in Hermitian mode
x = hifApply(hdl, b, 'SH');  % Hermitian triangular solve
b2 = hifApply(hdl, x, 'MH'); % Hermitian mat-vec multiplication
res = norm(b2 - b) / norm(b);
assert(res <= 1e-10, 'residual %g is too large', res);

delete(hdl);
