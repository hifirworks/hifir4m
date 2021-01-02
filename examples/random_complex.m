%{
This example uses a density 0.5 of 10x10 sparse matrix to show how to use
HILUCSI4M and the FGMRES solver.
%}
clear;
load('kim1.mat', 'Problem');
is_mixed = false;
is_complex = true;
A = Problem.A;
n = size(A,1);
b = A*ones(n,1);

%% Initialize database
h = HILUCSI(is_mixed, is_complex);

%% Factorize A
[tt, info] = h.factorize(A);
disp(info);
% test solve with 2 RHS
X = h.m_solve2([b b]);
assert(norm(X(:,1)-X(:,2))<=1e-12);

%% Solve for x=A\b
[x, flag] = h.fgmres(A, b);
if flag; fprintf(2, 'warning! solver failed with flag %d\n', flag); end
res = norm(A*x-b)/norm(b);
if res > 1e-6; fprintf(2, 'residual too large %.4g\n', res); end

%% Finalize
clear h