function rt = hifir4m_root
% get the root of the project

persistent root__
if isempty(root__); root__ = which('hifir4m_root'); end
rt = fileparts(root__);
end
