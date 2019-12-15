function rev = hilucsi4m_git_revision(pkg, rev_type)
%HILUCSI4M_GIT_REVISION - Get git revision of the package
%
% Syntax:
%   rev = hilucsi4m_git_revision
%   rev = hilucsi4m_git_revision(pkg)
%   rev = hilucsi4m_git_revision(___, rev_type)
%
% Description:
%   rev = hilucsi4m_git_revision provides the git revision of hilucsi4m
%   packages.
%
%   rev = hilucsi4m_git_revision(pkg) provides the git revision by given a
%   package, pkg can be either
%       'hilucsi4m' - matlab binding (default)
%       'hilucsi'   - C++ HILUCSI that the matlab binding was built upon
%
%   rev = hilucsi4m_git_revision(___, rev_type) provides options for
%   revision string types, rev_type must be either
%       'default' - default one (default)
%       'short'   - short revision hash
%
% Examples:
%   Get the default revision hash of HILUCSI4M package
%       >> disp(hilucsi4m_git_revision);  % should be a long hash
%
%   Get the short revision
%       >> disp(hilucsi4m_git_revision([], 'short'));
%
%   Get HILUCSI revision
%       >> disp(hilucsi4m_git_revision('hilucsi'));
%
% Note:
%   In case the command `git` fails, the returned revision string will be
%   'UNKNOWN'.
%
% See Also:
%   HILUCSI4M_ROOT

%------------------------- BEGIN MAIN CODE ------------------------------%

if nargin < 1 || isempty(pkg); pkg = 'hilucsi4m'; end
assert(ismember(pkg, {'hilucsi4m', 'hilucsi'}));
if nargin < 2 || isempty(rev_type); rev_type = 'default'; end
assert(ismember(rev_type, {'default', 'short'}));
cmd = 'git rev-parse HEAD';
if strcmp(rev_type, 'short'); cmd = 'git rev-parse --short HEAD'; end
d_bak = pwd;
cd(hilucsi4m_root);
if strcmp(pkg, 'hilucsi'); cd('hilucsi'); end
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