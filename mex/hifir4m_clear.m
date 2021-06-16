function hifir4m_clear(dbase)
%HIFIR4M_CLEAR - Clear the storage in a database
%
% Syntax:
%   hifir4m_clear(dbase)
%
% Description:
%   Clearing the storage is the purpose of this routine.
%
% See Also:
%   HIFIR4M_FINALIZE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: AGPLv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

hifir4m_mex(HIFIR4M_CLEAR, dbase);

%-------------------------- END MAIN CODE -------------------------------%
end