function hif = hifUpdate(hif, A)
% hifUpdate Updates a coefficient matrix in hif for iterative refinement
%
%   hif = hifUpdate(hif, A)
%
% See also hifApply hifRefactorize

if issparse(A)
    Astruct = hifir4m_sp2crs(A);
else
    Astruct = A;
    if hifir4m_isint64
        if ~isa(A.row_ptr, 'int64')
            Astruct.row_ptr = int64(Astruct.row_ptr);
            Astruct.col_ind = int64(Astruct.col_ind);
        end
    else
        if ~isa(A.row_ptr, 'int32')
            Astruct.row_ptr = int32(Astruct.row_ptr);
            Astruct.col_ind = int32(Astruct.col_ind);
        end
    end
end
if ~isempty(hif.A)
    assert(hif.A.nrows == Astruct.nrows, 'mismatched row sizes');
    assert(hif.A.ncols == Astruct.ncols, 'mismatched column sizes');
end
hif.A = Astruct;
end