function [y, varargout] = hifApply(hif, x, op, rank, nirs)
% hifApply  Applies a HIFIR preconditioner to solve or multiply.
%
%     y = hifApply(hif, x [op, rank, nirs])
% computes y = M^g x, y = M^gH x, y = M x, and y = M^H x for
% op = 'S' (default), 'SH' (or 'ST'), 'M', and 'MH' (or 'MT'),
% respectively. Rank specifies the rank of the final Schur
% complement; use 0 for the rank determined by hifCreate and
% use -1 for full rank (to machine precision). Nirs specifies
% the number of iterative refinements.
%
% See also hifCreate, Hifir

if nargin < 3
    op = 'S';
end
if nargin < 4
    rank = 0;
end
if nargin < 5
    nirs = 1;
end

if op(1) == 'S' || op(1) == 's'
    assert(nargin <= 5, 'Solve accepts up to five arguments.');
    if nargin == 5 && nirs > 1
        % with iterative refinement
        [y, varargout{:}] = hifir4m_mex(HifEnum.M_SOLVE, hif.hdl, ...
            A.row_ptr, A.col_ind, A.val, x, nirs, ...
            op(end) == 'H' || op(end) == 'T' || op(end) == 'h' || op(end) == 't', rank);
    else
        [y, varargout{1:nargout}] = hifir4m_mex(HifEnum.M_SOLVE, hif.hdl, x, ...
            op(end) == 'H' || op(end) == 'T' || op(end) == 'h' || op(end) == 't', rank);
    end
else
    assert(op(1) == 'M' || op(1) == 'm',  ...
        'First character must be ''s'' or ''m''');
    assert(nargin <= 4, 'Multiply accepts up to four arguments.');

    [y, varargout{1:nargout}] = hifir4m_mex(HifEnum.M_MULTIPLY, hif.hdl, x, ...
        op(end) == 'H' || op(end) == 'T' || op(end) == 'h' || op(end) == 't', rank);
end

end
