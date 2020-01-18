%{
Example of solving singular system
%}

clear;
h = HILUCSI;
N = 256;
A = gallery('neumann', 256);
left_nsp = null(full(A'));
b = rand(N, 1);
% filter the left null space from b, need to do it for now
b = b - ((left_nsp'*b)/norm(left_nsp)^2)*left_nsp;
h.factorize(A);
x = h.fgmres(A, b, [], [], [], [], [], [], [1, N]);
fprintf(1, 'relative residual in 2 norm is %g\n', norm(A*x-b)/norm(b));
clear h