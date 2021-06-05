function A = hifir4m_ensure_int(A)
%HIFIR4M_ENSURE_INT - Convert CRS to proper integer type (i.e., int{32,64})
%
% Syntax:
%   A = hifir4m_ensure_int(A)

isint64 = hifir4m_isint64;
if isa(A.row_ptr, 'int32') && isint64
    A.row_ptr = int64(A.row_ptr);
    A.col_ind = int64(A.col_ind);
elseif isa(A.row_ptr, 'int64') && ~isint64
    A.row_ptr = int32(A.row_ptr);
    A.col_ind = int32(A.col_ind);
end

end