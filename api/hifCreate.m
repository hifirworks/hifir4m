function [hdl, info, time] = hifCreate(A, S, varargin)
% hifCreate  Creates a HIFIR preconditioner for a given matrix
%
%    hdl = hifCreate(A [, S])
% creates a HIFIR preconditioner for a matrix A from A itself or from S,
% a sparsifier of A, and A and S can be MATLAB's sparse matrix or a CRS
% struct containing row_ptr, col_ind, vals entries. If S is empty, A will
% be used. The CRS structs are compatible with MATLAB Coder.
%
%    hdl = hifCreate(A, S, params)
% allows specifying the parameters using a HifParams struct. This mode is
% compatible with MATLAB Coder.
%
%    hdl = hifCreate(A, S, 'name1', val, 'name2', val2, ...)
% allows specifying the parameters using name-value pairs as for HifParams.
% This mode is not compatible with MATLAB Coder.
%
%    [hdl, info, time] = hifCreate(...)
% returns additional factorization and timing information.
%
% See also hifApply, Hifir, HifParams

% Copyright 2019--2021 NumGeom Group of Stony Brook University

%% Parse input arguments

if ~isempty(varargin) && isstruct(varargin{1})
    params = varargin{1};
else
    params = HifParams(varargin{:});
end

%% Create and setup the preconditioner
hif = hifir4m_mex(HifEnum.CREATE, params.is_mixed, params.is_complex);

% Setup object
if issparse(A)
    Astruct = hifir4m_sp2crs(A);
else
    if hifir4m_isint64
        if ~isa(A.row_ptr, 'int64')
            Astruct.row_ptr = int64(A.row_ptr);
            Astruct.col_ind = int64(A.col_ind);
        end
    elseif ~isint64
        if ~isa(A.row_ptr, 'int32')
            Astruct.row_ptr = int32(A.row_ptr);
            Astruct.col_ind = int32(A.col_ind);
        end
    end
    Astruct.val = double(A.val);
end

if nargin <= 1 || isempty(S)
    Sstruct = Astruct;
else
    if hifir4m_isint64
        if ~isa(S.row_ptr, 'int64')
            Sstruct.row_ptr = int64(S.row_ptr);
            Sstruct.col_ind = int64(S.col_ind);
        end
    elseif ~isint64
        if ~isa(S.row_ptr, 'int32')
            Sstruct.row_ptr = int32(S.row_ptr);
            Sstruct.col_ind = int32(S.col_ind);
        end
    end
    Sstruct.val = double(S.val);
end

[time, info] = hifir4m_mex(HifEnum.FACTORIZE, hif, ...
    Sstruct.row_ptr, Sstruct.col_ind, Sstruct.val, params);

%% Create HIFIR handle and return back to MATLAB
hdl = Hifir(hif, params, Astruct);

end
