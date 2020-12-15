function A = hilucsi4m_2int64(A)
%HILUCSI4M_2INT64 - Convert int32 CRS to int64 CRS
%
% Syntax:
%   A = hilucsi4m_2int64(A)
A.row_ptr = int64(A.row_ptr);
A.col_ind = int64(A.col_ind);
end