function b = crs_prodAx(A, x, b)
%crs_prodAx Compute b=A*x for a sparse matrix A in CRS format,
% assuming that number of columns of A is equal to size(x,1).
%
%      b = crs_prodAx(A, x [,b])
% Computes b=A*x in serial.
%
% See http://www.netlib.org/linalg/html_templates/node98.html

if nargin == 2
    b = coder.nullcopy(zeros(A.nrows, size(x, 2)));
else
    if size(b, 1) < A.nrows || size(b, 2) < size(x, 2)
        error('crs_prodAx:BufferTooSmal', 'Buffer space for output b is too small.');
    end
end

%% Compute b=A*x
b = crs_prodAx_kernel(A.row_ptr, A.col_ind, A.val, x, int32(size(x, 1)), ...
    b, int32(size(b, 1)), A.nrows, int32(size(x, 2)));

function b = crs_prodAx_kernel(row_ptr, col_ind, val, ...
    x, x_m, b, b_m, nrows, nrhs)

istart = int32(1); iend = nrows;

xoffset = int32(0); boffset = int32(0);
for k = 1:nrhs
    for i = istart:iend
        t = 0.0;
        for j = row_ptr(i):row_ptr(i+1) - 1
            t = t + val(j) * x(xoffset+col_ind(j));
        end
        b(boffset+i) = t;
    end
    xoffset = xoffset + x_m; boffset = boffset + b_m;
end

