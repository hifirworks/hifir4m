function [rowptr, colind, vals] = hilucsi4m_sp2crs(A)
%HILUCSI4M_SP2CRS - Convert sparse to CRS representation

assert(issparse(A), 'input A must be sparse');
[m, n] = size(A);
if m ~= n; error('input matrix must be squared!'); end
[rs, cs, vs] = find(A);
rs = int32(rs);
cs = int32(cs);
[rowptr, colind, vals] = hilucsi4m_ijv2crs(int32(m), rs, cs, vs);
end
