function hilucsi4m_finalize(dbase)
%HILUCSI4M_FINALIZE - Finalize a low-level database
%
% Syntax:
%   hilucsi4m_finalize(dbase)
%
% Description:
%   Clean up an existing database.
%
% See Also:
%   HILUCSI4M_INITIALIZE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: GLPv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

hilucsi4m_mex(HILUCSI4M_DESTROY, dbase);

%-------------------------- END MAIN CODE -------------------------------%
end