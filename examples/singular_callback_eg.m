%{
Example of solving singular system

NOTE:
    This example uses user callback to filter out the constant mode (right
    null space). The callback must be taking one argument and return one
    argument, where the argument is a 1D column vector, which is the
    direction vector in GMRES. The result is expected to not contain any
    components in the (right) null space.
%}

clear;
h = HIF;
N = 256;
A = gallery('neumann', 256);
left_nsp = null(full(A'));
b = rand(N, 1);
% filter the left null space from b, need to do it for now
b = b - ((left_nsp'*b)/norm(left_nsp)^2)*left_nsp;
h.factorize(A);
% one can do this
x = h.gmres(A, b, [], [], [], [], [], [], @my_filter);
% or
% x = h.fgmres(A, b, [], [], [], [], [], [], @(x) my_filter(x, N));
fprintf(1, 'relative residual in 2 norm is %g\n', norm(A*x-b)/norm(b));
clear h

function v = my_filter(v, N)
% NOTE and WARING: The function must agree with this interface, input is
% one column vector with output one column vector. This will be called
% inside the GMRES KSP solver. Notice that one can capture user data by
% having multiple arguments where default values are supplied from the
% second one, or create function handle by capturing other data vars.
if nargin < 2; N = length(v); end
v = v - sum(v)/N;
end