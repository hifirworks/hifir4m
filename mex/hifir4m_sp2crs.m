function A_struct = hifir4m_sp2crs(A)
%HIFIR4M_SP2CRS - Convert sparse to CRS representation
%
% Syntax:
%   A = hifir4m_sp2crs(A)
%
% Description:
%   Convert sparse matrix into a struct of CRS with `row_ptr`, `col_ind`
%   and `val`

assert(issparse(A), 'input A must be sparse');
[m, n] = size(A);
if m ~= n; error('input matrix must be squared!'); end
[rs, cs, vs] = find(A);
isint64 = hifir4m_isint64;
if ~isint64
    rs = int32(rs);
    cs = int32(cs);
else
    rs = int64(rs);
    cs = int64(cs);
end
[A_struct.row_ptr, A_struct.col_ind, A_struct.val] = hifir4m_ijv2crs(m, rs, ...
    cs, vs);
end
