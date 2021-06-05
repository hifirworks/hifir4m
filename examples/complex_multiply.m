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
[tt, info] = h.factorize(A);
disp(info);

%% Solve for x=h\b
x = h.solve(b);

%% Compute y = h*x
y = h.mmultiply(x);

assert(norm(y-b)/norm(b) <= 1e-12, 'failed matrix-vector test');

%% Solve for x=h'\b
x = h.solve(b, true);

%% Compute y = h'*x
y = h.mmultiply(x, true);

assert(norm(y-b)/norm(b) <= 1e-12, 'failed Hermitian matrix-vector test');

%% Finalize
clear h