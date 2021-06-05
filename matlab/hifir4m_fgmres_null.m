function varargout = hifir4m_fgmres_null(dbase, A, b, varargin)
%HIFIR4M_FGMRES_NULL - right-preconditioned FGMRES for null-space vectors
%
% Syntax:
%   x = hifir4m_fgmres_null(dbase, A, b)
%
% Long description
%
% See Also:
%   HIFIR_FGMRES

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: GLPv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

gmres_pars = [30 1e-13 500]; % restart,rtol,maxit
x0 = [];
verbose = true;
hiprec = false;
for i = 1:min(3, length(varargin))
    if ~isempty(varargin{i}); gmres_pars(i) = varargin{i}; end
end
if length(varargin) > 3; x0 = varargin{4}; end
if length(varargin) > 4
    if ~isempty(varargin{5}); verbose = logical(varargin{5}); end
end
if length(varargin) > 5
    if ~isempty(varargin{6}); hiprec = logical(varargin{6}); end
end
if isempty(x0)
    if isreal(b)
        x0 = zeros(size(b));
    else
        x0 = complex(zeros(size(b)), zeros(size(b)));
    end
end
% Convert to zero-based CRS
if issparse(A); A = hifir4m_sp2crs(A); end
A = hifir4m_ensure_int(A);
[varargout{1:nargout}] = hifir4m_mex(HIFIR4M_KSP_NULL, dbase, ...
    A.row_ptr, A.col_ind, A.val, b, gmres_pars(1), gmres_pars(2), ...
    gmres_pars(3), x0, verbose, hiprec);

%-------------------------- END MAIN CODE -------------------------------%
end