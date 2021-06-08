function opts = hifir4m_create_params(varargin)
%HIFIR4M_CREATE_PARAMS - Create options for HIFIR preconditioner
%
% Syntax:
%   opts = hifir4m_create_params
%   opts = hifir4m_create_params(OptionFieldName, OptionValue)
%
% Description:
%   HIFIR4M_CREATE_PARAMS constructs the control parameters used in
%   factorizing HIFIR preconditioner. The parameters are well-documented
%   in the c++ code thus omitting here.
%
%   opts = hifir4m_create_params gets the default parameters
%
%   opts = hifir4m_create_params(OptionFieldName, OptionValue) mimics
%   the MATLAB parameter syntax (key, value), which is commonly used (e.g.
%   in 2-D line plotting). Available option field names are
%       tau_L, tau_U, tau_d, tau_kappa, alpha_L, alpha_U, rho, c_d,
%       c_h, N, verbose, rf_par, reorder, saddle, check, pre_scale,
%       symm_pre_lvls, threads,
%   whose meanings are documented in the c++ code.
%
% Examples:
%   To get default constrol parameters
%       >> opts = hifir4m_create_params;
%
%   To customize parameters
%       >> opts = hifir4m_create_params('tau_L', 1e-3, 'tau_U', 1e-3);
%       >> assert(opts.tau_L == 1e-3);
%
% See Also:
%   HIFIR4M_FACTORIZE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: GLPv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

persistent fnames
if isempty(fnames)
    fnames = {'tau_L', 'tau_U', 'kappa_d', 'kappa', 'alpha_L', 'alpha_U', ...
        'rho', 'c_d', 'c_h', 'N', 'verbose', 'rf_par', 'reorder', 'saddle', ...
        'check', 'pre_scale', 'symm_pre_lvls', 'threads', 'mumps_blr', ...
        'fat_schur_1st', 'rrqr_cond', 'pivot', 'gamma', 'beta', 'is_symm', ...
        'no_pre'};
end

p = inputParser;
addOptional(p, 'tau_L', 1e-4, @(x) isscalar(x) && x > 0);
addOptional(p, 'tau_U', 1e-4, @(x) isscalar(x) && x > 0);
addOptional(p, 'kappa_d', 3.0, @(x) isscalar(x) && x > 0);
addOptional(p, 'kappa', 3.0, @(x) isscalar(x) && x > 0);
addOptional(p, 'alpha_L', 10, @(x) isscalar(x) && x > 0);
addOptional(p, 'alpha_U', 10, @(x) isscalar(x) && x > 0);
addOptional(p, 'rho', 0.5, @(x) isscalar(x) && x > 0);
addOptional(p, 'c_d', 10.0, @(x) isscalar(x) && x > 0);
addOptional(p, 'c_h', 2.0, @(x) isscalar(x) && x > 0);
addOptional(p, 'N', -1, @isscalar);
addOptional(p, 'verbose', 1, @(x) isscalar(x) && x >= 0);
addOptional(p, 'rf_par', 1, @(x) isscalar(x) && x >= 0);
addOptional(p, 'reorder', 1, @(x) isscalar(x) && x >= 0);
addOptional(p, 'saddle', 1, @(x) isscalar(x) && x >= 0);
addOptional(p, 'check', 1, @(x) isscalar(x) && x >= 0);
addOptional(p, 'pre_scale', 0, @(x) isscalar(x) && x >= 0);
addOptional(p, 'symm_pre_lvls', 1, @(x) isscalar(x) && x >= 0);
addOptional(p, 'threads', 0, @(x) isscalar(x) && x >= 0);
addOptional(p, 'mumps_blr', 1, @(x) isscalar(x) && x >= 0);
addOptional(p, 'fat_schur_1st', 0, @(x) isscalar(x) && x >= 0);
addOptional(p, 'rrqr_cond', 0, @(x) isscalar(x));
addOptional(p, 'pivot', 0, @(x) isscalar(x) && x >= 0);
addOptional(p, 'gamma', 1.0, @(x) isscalar(x) && x >= 0);
addOptional(p, 'beta', 1e3, @(x) isscalar(x) && x >= 0);
addOptional(p, 'is_symm', 0, @(x) isscalar(x) && x>= 0);
addOptional(p, 'no_pre', 0, @(x) isscalar(x) && x >= 0);
parse(p, varargin{:});
sorted_opts = p.Results;
% NOTE parser gives sorted structure
opts = struct;
for idx = 1:length(fnames)
    f = fnames{idx};
    opts = setfield(opts, f, getfield(sorted_opts, f));
end

%-------------------------- END MAIN CODE -------------------------------%
end
