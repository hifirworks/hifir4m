clear;
load('kim1.mat', 'Problem');
is_mixed = false;
is_complex = true;
A = Problem.A;
n = size(A,1);
b = A*ones(n,1);

%% Initialize database
h = HIF(is_mixed, is_complex);

%% Factorize A
info = h.factorize(A);
disp(info);

%% Solve for x=A\b
[x, flag] = h.fgmres(A, b);
if flag; fprintf(2, 'warning! solver failed with flag %d\n', flag); end
res = norm(A*x-b)/norm(b);
if res > 1e-6; fprintf(2, 'residual too large %.4g\n', res); end

%% Finalize
clear h