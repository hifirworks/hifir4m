function A_struct = hilucsi4m_sp2crs(A)
%HILUCSI4M_SP2CRS - Convert sparse to CRS representation
%
% Syntax:
%   A = hilucsi4m_sp2crs(A)
%
% Description:
%   Convert sparse matrix into a struct of CRS with `row_ptr`, `col_ind`
%   and `val`

assert(issparse(A), 'input A must be sparse');
[m, n] = size(A);
if m ~= n; error('input matrix must be squared!'); end
[rs, cs, vs] = find(A);
rs = int32(rs);
cs = int32(cs);
[A_struct.row_ptr, A_struct.col_ind, A_struct.val] = hilucsi4m_ijv2crs(...
    int32(m), rs, cs, vs);
end
