function rt = hilucsi4m_root
% get the root of the project

persistent root__
if isempty(root__); root__ = which('hilucsi4m_root'); end
rt = fileparts(root__);
end