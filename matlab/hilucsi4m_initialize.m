function dbase = hilucsi4m_initialize(is_mixed)
%HILUCSI4M_INITIALIZE - Initialize a low-level database
%
% Syntax:
%   dbase = hilucsi4m_initialize
%   dbase = hilucsi4m_initialize(true)
%
% Description:
%   HILUCSI4M_INITIALIZE should be the first routine you call in this
%   library, as it performs the initialization of an internal database.
%
%   dbase = hilucsi4m_initialize constructs a default instance of database
%   with double precision preconditioner.
%
%   dbase = hilucsi4m_initialize(true) constructs a mixed-precision
%   instance of the internal database, i.e., the preconditioner uses single
%   precision floating numbers whereas the solvers accept double precision.
%
% See Also:
%   HILUCSI4M_FINALIZE, HILUCSI4M_FACTORIZE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: GLPv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

if nargin < 1 || isempty(is_mixed); is_mixed = false; end
assert(isscalar(is_mixed), 'is_mixed must be scalar (prefer boolean)');
dbase = hilucsi4m_mex(HILUCSI4M_CREATE, is_mixed);

%-------------------------- END MAIN CODE -------------------------------%
end