function rev = hifir4m_git_revision(pkg, rev_type)
%HIFIR4M_GIT_REVISION - Get git revision of the package
%
% Syntax:
%   rev = hifir4m_git_revision
%   rev = hifir4m_git_revision(pkg)
%   rev = hifir4m_git_revision(___, rev_type)
%
% Description:
%   rev = hifir4m_git_revision provides the git revision of hifir4m
%   packages.
%
%   rev = hifir4m_git_revision(pkg) provides the git revision by given a
%   package, pkg can be either
%       'hifir4m' - matlab binding (default)
%       'hifir'   - C++ HIFIR that the matlab binding was built upon
%
%   rev = hifir4m_git_revision(___, rev_type) provides options for
%   revision string types, rev_type must be either
%       'default' - default one (default)
%       'short'   - short revision hash
%
% Examples:
%   Get the default revision hash of HIFIR4M package
%       >> disp(hifir4m_git_revision);  % should be a long hash
%
%   Get the short revision
%       >> disp(hifir4m_git_revision([], 'short'));
%
%   Get HIFIR revision
%       >> disp(hifir4m_git_revision('hifir'));
%
% Note:
%   In case the command `git` fails, the returned revision string will be
%   'UNKNOWN'.
%
% See Also:
%   HIFIR4M_ROOT

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: AGPLv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

if nargin < 1 || isempty(pkg); pkg = 'hifir4m'; end
assert(ismember(pkg, {'hifir4m', 'hifir'}));
if nargin < 2 || isempty(rev_type); rev_type = 'default'; end
assert(ismember(rev_type, {'default', 'short'}));
cmd = 'git rev-parse HEAD';
if strcmp(rev_type, 'short'); cmd = 'git rev-parse --short HEAD'; end
d_bak = pwd;
cd(hifir4m_root);
if strcmp(pkg, 'hifir'); cd('hifir'); end
[flag, rev] = system(cmd);
if flag
    fprintf(2, 'Could not get git revision (%d)!\n', flag);
    rev = 'UNKNOWN';
else
    rev = strip(rev);
end
cd(d_bak);

%-------------------------- END MAIN CODE -------------------------------%
end