function dbase = hifir4m_initialize(is_mixed, is_complex)
%HIFIR4M_INITIALIZE - Initialize a low-level database
%
% Syntax:
%   dbase = hifir4m_initialize
%   dbase = hifir4m_initialize(true)
%   dbase = hifir4m_initialize(___, true)
%
% Description:
%   HIFIR4M_INITIALIZE should be the first routine you call in this
%   library, as it performs the initialization of an internal database.
%
%   dbase = hifir4m_initialize constructs a default instance of database
%   with double precision preconditioner.
%
%   dbase = hifir4m_initialize(true) constructs a mixed-precision
%   instance of the internal database, i.e., the preconditioner uses single
%   precision floating numbers whereas the solvers accept double precision.
%
%   dbase = hifir4m_initialize(___, true) constructs a complex number
%   preconditioner instance.
%
% See Also:
%   HIFIR4M_FINALIZE, HIFIR4M_FACTORIZE

% Author: Qiao Chen
% Email: qiao.chen@stonybrook.edu
% License: GLPv3+

%------------------------- BEGIN MAIN CODE ------------------------------%

if nargin < 1 || isempty(is_mixed); is_mixed = false; end
if nargin < 2 || isempty(is_complex); is_complex = false; end
assert(isscalar(is_mixed), 'is_mixed must be scalar (prefer boolean)');
dbase = hifir4m_mex(HIFIR4M_CREATE, is_mixed, is_complex);

%-------------------------- END MAIN CODE -------------------------------%
end