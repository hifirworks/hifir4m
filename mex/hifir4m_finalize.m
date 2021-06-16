function hifir4m_finalize(dbase)
%HIFIR4M_FINALIZE - Finalize a low-level database
%
% Syntax:
%   hifir4m_finalize(dbase)
%
% Description:
%   Clean up an existing database.
%
% See Also:
%   HIFIR4M_INITIALIZE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: AGPLv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

hifir4m_mex(HIFIR4M_DESTROY, dbase);

%-------------------------- END MAIN CODE -------------------------------%
end