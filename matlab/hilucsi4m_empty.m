function flag = hilucsi4m_empty(dbase)
%HILUCSI4M_EMPTY - Check the emptyness of a database
%
% Syntax:
%   hilucsi4m_empty(dbase)
%
% Description:
%   Check the emptyness of a database
%
% See Also:
%   HILUCSI4M_INITIALIZE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: GLPv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

flag = hilucsi4m_mex(HILUCSI4M_CHECK, dbase);

%-------------------------- END MAIN CODE -------------------------------%
end