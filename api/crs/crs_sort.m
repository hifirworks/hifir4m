function [col_ind, val] = crs_sort(row_ptr, col_ind, val)
%crs_sort  Sort column indices within each row..
%  [col_ind, val] = crs_sort(row_ptr, col_ind, val)

%#codegen -args {coder.typeof(int32(0), [inf,1]),
%#codegen coder.typeof(int32(0), [inf,1]), coder.typeof(0, [inf,1])}
%#codegen crs_sort0 -args {coder.typeof(int32(0), [inf,1]),
%#codegen coder.typeof(int32(0), [inf, 1])}

assert(nargin==nargout+1 && nargin>=2);

for i=1:int32(length(row_ptr))-1
    ascend = true;
    for j=row_ptr(i)+1 : row_ptr(i+1)-1
        if col_ind(j)<col_ind(j-1)
            ascend = false;
            break;
        end
    end
    
    if ~ascend
        if nargin==3
            % copy data into buffer
            buf_indx = zeros(row_ptr(i+1)-row_ptr(i),1,'int32');
            buf_val = zeros(row_ptr(i+1)-row_ptr(i),1, class(val));
            ind = 1;
            for j=row_ptr(i) : row_ptr(i+1)-1
                buf_indx(ind) = col_ind(j);
                buf_val(ind) = val(j);
                ind=ind+1;
            end
            
            % sort in place
            [buf_indx, buf_val] = heapsort_tag(buf_indx, buf_val);
            
            % copy data bacl fom buffer
            ind = 1;
            for j=row_ptr(i) : row_ptr(i+1)-1
                col_ind(j) = buf_indx(ind);
                val(j) = buf_val(ind);
                ind=ind+1;
            end
        else
            % copy data into buffer
            buf_indx = coder.nullcopy(zeros(row_ptr(i+1)-row_ptr(i),1,'int32'));
            ind = 1;
            for j=row_ptr(i) : row_ptr(i+1)-1
                buf_indx(ind) = col_ind(j);
                ind=ind+1;
            end
            
            % sort in place
            buf_indx = heapsort(buf_indx);
            
            % copy data bacl fom buffer
            ind = 1;
            for j=row_ptr(i) : row_ptr(i+1)-1
                col_ind(j) = buf_indx(ind);
                ind=ind+1;
            end
        end
    end
end
