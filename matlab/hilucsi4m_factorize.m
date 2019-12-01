function varargout = hilucsi4m_factorize(dbase, A, opts)
%HILUCSI4M_FACTORIZE - Factorize HILUCSI preconditioner

if nargin < 3 || isempty(opts); opts = hilucsi4m_create_options; end
% convert A to zero-based CRS
[rowptr, colind, vals] = hilucsi4m_sp2crs(A);
[varargout{1:nargout}] = hilucsi4m_mex(HILUCSI4M_FACTORIZE, dbase, rowptr, ...
    colind, vals, opts);
end