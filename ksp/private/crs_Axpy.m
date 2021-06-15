function y = crs_Axpy(A, x, y, nthreads)
%crs_Axpy Compute y=y+A*x for a sparse matrix A in CRS format,
% assuming that number of columns of A is equal to size(x,1).
%
%      y = crs_Axpy(A, x ,y)
% Computes y=A*x+y in serial.
%
%      y = crs_Axpy(A, x, y, nthreads)
% Computes y=A*x+y Testing nthreads OpenMP threads, where nthreads is an int32.
% In this mode, there are implicit barriers in the function.
%
%      y = crs_Axpy(A, x, y, [])
% Computes y=A*x+y locally within each OpenMP thread, assuming the
% parallel region has already been initialized. The argument y has
% been preallocated. In this mode, there is no implicit barrier.
%
%      y = crs_Axpy(A, x, y, nthreads, comm)
% also performs communication within the MPI communicator.
% It assumes that A is partitioned columnwise, and hence A*x is
% a partial sum on each process. An allreduce is performed on y
% within the communicator.
%
% Note that if A is partitioned rowwise, then y must be duplicated
% on all processes, and prodAx should be called without the communicator.
%
%      y = crs_Axpy(A, x, y, nthreads, comm, pgmsg [, pgsz])
% specifies an additional message that should be summed during
% allreduce and strored in y. This is to piggy-back on allreduce
% to help aggregate small messages to optimize MPI communication.
%
% The compiled version reqiures A and x to be double precision.
% The M file works for both single and double precision.

% See http://www.netlib.org/linalg/html_templates/node98.html

%#codegen -args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),int32(0)}
%#codegen crs_Axpy_ser1 -args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}

if size(y, 1) < A.nrows || size(y, 2) < size(x, 2)
    m2c_error('crs_Axpy:BufferTooSmal', 'Buffer space for output y is too small.');
end


if nargin > 3 && ~isempty(nthreads)
    %% Declare parallel region
    if ~ompGetNested && ompGetNumThreads > 1 && nthreads(1) > 1
        OMP_begin_master
        m2c_warn('crs_Axpy:NestedParallel', ...
            'You are trying to use nested parallel regions. Solution may be incorrect.');
        OMP_end_master
    end

    OMP_begin_parallel(int32(nthreads(1)));

    %% Compute y=A*x+y
    y = crs_Axpy_kernel(A.row_ptr, A.col_ind, A.val, x, int32(size(x, 1)), ...
        y, int32(size(y, 1)), A.nrows, int32(size(x, 2)), ompGetNumThreads > 1);

    %% End parallel region
    OMP_end_parallel;
elseif nargin < 4 || ompGetNumThreads == 1
    %% Compute y=A*x+y
    y = crs_Axpy_kernel(A.row_ptr, A.col_ind, A.val, x, int32(size(x, 1)), ...
        y, int32(size(y, 1)), A.nrows, int32(size(x, 2)));
else
    %% Compute y=A*x+y
    y = crs_Axpy_kernel(A.row_ptr, A.col_ind, A.val, x, int32(size(x, 1)), ...
        y, int32(size(y, 1)), A.nrows, int32(size(x, 2)), true);
end

function y = crs_Axpy_kernel(row_ptr, col_ind, val, ...
    x, x_m, y, b_m, nrows, nrhs, varargin)

coder.inline('never');
if ~isempty(varargin) && varargin{1}
    [istart, iend] = OMP_local_chunk(nrows);
else
    istart = int32(1); iend = nrows;
end

xoffset = int32(0); boffset = int32(0);
for k = 1:nrhs
    for i = istart:iend
        t = y(boffset+i);

        for j = row_ptr(i):row_ptr(i+1) - 1
            t = t + val(j) * x(xoffset+col_ind(j));
        end
        y(boffset+i) = t;
    end
    xoffset = xoffset + x_m; boffset = boffset + b_m;
end

function test %#ok<DEFNU>
%!test
%! if ~exist(['crs_Axpy.' mexext], 'file')
%!    m=100; n = 20;
%! else
%!    m=10000; n = 2000;
%! end
%! tic; sp = sprand(m,n,0.5); x=rand(size(sp,2),2);
%! fprintf(1, '\n\tGenerated random matrix in %g seconds\n', toc);
%! tic; b0 = sp*x + ones(m, 2);
%! fprintf(1, '\tComputed reference solution in %g seconds\n', toc);
%! tic; A = crs_matrix(sp);
%! fprintf(1, '\tConverted into crs_matrix in %g seconds\n', toc);
%! fprintf(1, '\tTesting serial: ');
%! tic; b1 = crs_Axpy(A, x, ones(m, 2));
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert(norm(b0-b1)<=1.e-10);

%! b2 = zeros(size(sp,1),2);
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>ompGetMaxThreads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     tic; b2 = crs_Axpy(A, x, ones(m, 2), nthreads);
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert(norm(b0-b2)<=1.e-10);
%! end
