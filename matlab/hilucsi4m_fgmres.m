function varargout = hilucsi4m_fgmres(dbase, A, b, varargin)
%HILUCSI4M_FGMRES - Flexible (right-preconditioned) GMRES with HILUCSI

gmres_pars = [30 1e-6 500]; % restart,rtol,maxit
x0 = [];
verbose = true;
for i = 1:min(3, length(varargin)); gmres_pars(i) = varargin{i}; end
if length(varargin) > 3; x0 = varargin{4}; end
if length(varargin) > 4; verbose = logical(varargin{5}); end
if isempty(x0); x0 = zeros(size(b)); end
% Convert to zero-based CRS
[rowptr, colind, vals] = hilucsi4m_sp2crs(A);
[varargout{1:nargout}] = hilucsi4m_mex(HILUCSI4M_KSP_SOLVE, dbase, rowptr, ...
    colind, vals, b, gmres_pars(1), gmres_pars(2), gmres_pars(3), x0, verbose);
end