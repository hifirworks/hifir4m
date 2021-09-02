function varargout = hifRefactorize(hif, S, varargin)
% hifRefactorize refactorizes a preconditioner given a sparsifier S
%
%   hifRefactorize(hif, S)
%   hifRefactorize(hif, S, 'alpha_L', 5, 'alpha_U', 5);
%
% See also hifUpdate hifCreate hifApply

assert(~isempty(hif.hdl), 'the preconditioner cannot be empty.');

if ~isempty(varargin) && isstruct(varargin{1})
    params = varargin{1};
else
    params = HifParams(varargin{:});
end

% Setup sparsifier
if issparse(S)
    Sstruct = hifir4m_sp2crs(S);
else
    Sstruct = S;
    if hifir4m_isint64
        if ~isa(S.row_ptr, 'int64')
            Sstruct.row_ptr = int64(S.row_ptr);
            Sstruct.col_ind = int64(S.col_ind);
        end
    else
        if ~isa(S.row_ptr, 'int32')
            Sstruct.row_ptr = int32(S.row_ptr);
            Sstruct.col_ind = int32(S.col_ind);
        end
    end
end

[varargout{1:nargout-1}] = hifir4m_mex(HifEnum.FACTORIZE, hif.hdl, ...
    Sstruct.row_ptr, Sstruct.col_ind, Sstruct.val, params);

end