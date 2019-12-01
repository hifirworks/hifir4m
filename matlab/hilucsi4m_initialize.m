function dbase = hilucsi4m_initialize(is_mixed)
%HILUCSI4M_INITIALIZE - Initialize a low-level database

if nargin < 1 || isempty(is_mixed); is_mixed = false; end
assert(isscalar(is_mixed), 'is_mixed must be scalar (prefer boolean)');
dbase = hilucsi4m_mex(HILUCSI4M_CREATE, is_mixed);
end