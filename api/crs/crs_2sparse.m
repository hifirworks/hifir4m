function sp = crs_2sparse(varargin)
%crs_2sparse Create a MATLAB sparse matrix from CRS format.
%      sp = crs_2sparse(row_ptr, col_ind, val [,m ,n])
%      sp = crs_2sparse(A [,m ,n])
%
% See also crs_matrix
%
% Note: This function is incompatible with MATLAB Coder.
% It is provided for convenience of testing.
% It does not use multithreading.

assert(nargin==1 || nargin==3 || nargin==5);

if isstruct(varargin{1})
    row_ptr = varargin{1}.row_ptr;
    col_ind = varargin{1}.col_ind;
    val = varargin{1}.val;
    
    row_ind = crs_rowind(row_ptr, col_ind);
    
    if nargin==3
        sp = sparse(double(row_ind), double(col_ind), double(val), ...
            double(varargin{2}), double(varargin{3}));
    else
        sp = sparse(double(row_ind), double(col_ind), double(val), ...
            double(varargin{1}.nrows), double(varargin{1}.ncols));
    end
else
    row_ptr = varargin{1};
    col_ind = varargin{2};
    val = varargin{3};
    row_ind = crs_rowind(row_ptr, col_ind);
    
    if nargin==3
        sp = sparse(double(row_ind), double(col_ind), double(val));
    else
        sp = sparse(double(row_ind), double(col_ind), double(val), ...
            double(varargin{4}), double(varargin{5}));
    end
end

%!test
%! for k=1:100
%!     sp = sprand(20,10,0.5); x = rand(10,1);
%!     A = crs_matrix(sp);
%!     sp1 = crs_2sparse(A, int32(20), int32(10));
%!
%!     [row_ptr, col_ind, val] = crs_matrix(sp);
%!     sp1 = crs_2sparse(row_ptr, col_ind, val, int32(20), int32(10));
%!
%!     assert(isequal(sp, sp1));
%! end
