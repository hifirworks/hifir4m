function varargout = hilucsi4m_m_solve(dbase, b)
%HILUCSI4M_M_SOLVER - Accessing inv(M)

[varargout{1:nargout}] = hilucsi4m_mex(HILUCSI4M_M_SOLVE, dbase, b);
end