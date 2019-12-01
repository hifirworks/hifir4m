function hilucsi4m_finalize(dbase)
%HILUCSI4M_FINALIZE - Finalize a low-level database

hilucsi4m_mex(HILUCSI4M_DESTROY, dbase);
end