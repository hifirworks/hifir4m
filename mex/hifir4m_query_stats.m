function info = hifir4m_query_stats(dbase)
%HIFIR4M_QUERY_STATS - Query factorization statistics
%
% Syntax:
%   info = hifir4m_query_stats(dbase)
%
% Description:
%   Query factorization stats, such as nnz, levels, rank, etc. The output
%   info is a structure containing information fields.
%
% See Also:
%   HIFIR4M_FACTORIZE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: AGPLv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

info = hifir4m_mex(HIFIR4M_QUERY, dbase);

%-------------------------- END MAIN CODE -------------------------------%

end