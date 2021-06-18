%{
    Example of solving singular system using PIPIT with testing matrix
    gallery('neumann', 256).
%}

clear;

A = gallery('neumann', 256);
b = rand(size(A, 1), 1);

% Call PIPIT with 1 dim(Null)
x = pipitHifir(A, b, 1);
res = norm(A' * (b - A* x)) / norm(A' * b);
if res > 1e-6
    fprintf(2, 'residual too large %.4g\n', res);
end
