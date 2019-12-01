%{
This example uses a density 0.5 of 10x10 sparse matrix to show how to use
HILUCSI4M and the FGMRES solver.
%}
clear;

is_mixed = true;

A = sprand(10, 10, 0.5);
b = rand(10, 1);

%% Initialize database
dbase = hilucsi4m_initialize(is_mixed);

%% Factorize A
[tt, info] = hilucsi4m_factorize(dbase, A);
disp(info);

%% Solve for x=A\b
[x, flag] = hilucsi4m_fgmres(dbase, A, b);
if flag; fprintf(2, 'warning! solver failed with flag %d\n', flag); end
res = norm(A*x-b)/norm(b);
if res > 1e-6; fprintf(2, 'residual too large %.4g\n', res); end

%% Finalize
hilucsi4m_finalize(dbase);
clear dbase