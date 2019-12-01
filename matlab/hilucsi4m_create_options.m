function opts = hilucsi4m_create_options(varargin)
%HILUCSI4M_CREATE_OPTIONS - Create options for HILUCSI preconditioner

persistent fnames
if isempty(fnames)
    fnames = {'tau_L', 'tau_U', 'tau_d', 'tau_kappa', 'alpha_L', 'alpha_U', ...
        'rho', 'c_d', 'c_h', 'N', 'verbose', 'rf_par', 'reorder', 'saddle', ...
        'check', 'pre_scale', 'symm_pre_lvls', 'threads'};
end

p = inputParser;
addOptional(p, 'tau_L', 1e-4, @(x) isscalar(x) && x > 0);
addOptional(p, 'tau_U', 1e-4, @(x) isscalar(x) && x > 0);
addOptional(p, 'tau_d', 3.0, @(x) isscalar(x) && x > 0);
addOptional(p, 'tau_kappa', 3.0, @(x) isscalar(x) && x > 0);
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
parse(p, varargin{:});
sorted_opts = p.Results;
% NOTE parser gives sorted structure
opts = struct;
for idx = 1:length(fnames)
    f = fnames{idx};
    opts = setfield(opts, f, getfield(sorted_opts, f));
end
end