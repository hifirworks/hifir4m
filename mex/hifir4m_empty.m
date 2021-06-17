function flag = hifir4m_empty(dbase)
%HIFIR4M_EMPTY - Check the emptyness of a database
%
% Syntax:
%   hifir4m_empty(dbase)
%
% Description:
%   Check the emptyness of a database
%
% See Also:
%   HIFIR4M_INITIALIZE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: AGPLv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

flag = hifir4m_mex(HIFIR4M_CHECK, dbase);

%-------------------------- END MAIN CODE -------------------------------%
end