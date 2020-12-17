classdef HILUCSI
    %HILUCSI - OOP class representation of HILUCSI4M project
    %
    % Description:
    %   This class provides another way to use this package in an OOP
    %   fashion, which can be benefitial for certain programmers.

    properties (Access = protected)
        dbase  % internal database structure
    end

    methods
        function obj = HILUCSI(varargin)
            %HILUCSI - Construct an instance with internal database
            %
            % Syntax:
            %   hl = HILUCSI
            %   hl = HILUCSI(true)
            %   hl = HILUCSI(___, true)
            %
            % Description:
            %   hl = HILUCSI simpily constructs an instance without mixed
            %   precision support.
            %
            %   hl = HILUCSI(true) in constrast to the constructor above,
            %   this one enables mixed precision.
            %
            % See Also: HILUCSI4M_INITIALIZE
            obj.dbase = hilucsi4m_initialize(varargin{:});
        end
        function flag = is_mixed(obj)
            %IS_MIXED - Check if the HILUCSI instance is mixed precision
            %
            % See Also: ID
            flag = obj.dbase.is_mixed;
        end
        function flag = is_real(obj)
            %IS_REAL - Check if the HILUCSI instance is real arithmetic
            %
            % See Also: ID
            flag = obj.dbase.is_real;
        end
        function id_ = id(obj)
            %ID - Get the ID tag of the database
            %
            % See Also: IS_MIXED
            id_ = obj.dbase.id;
        end        
        function flag = empty(obj)
            %EMPTY - Check emptyness
            %
            % See Also: ID
            flag = hilucsi4m_empty(obj.dbase);
        end
        function varargout = factorize(obj, A, varargin)
            %FACTORIZE - Perform MLILU factorization
            %
            % See Also: HILUCSI4M_FACTORIZE, M_SOLVE
            [varargout{1:nargout}] = hilucsi4m_factorize(obj.dbase, A, ...
                varargin{:});
        end
        function hilu = export(obj, varargin)
            %EXPORT - Export internal data in C++ to MATLAB
            %
            % See Also: HILUCSI4M_EXPORT
            hilu = hilucsi4m_export(obj.dbase, varargin{:});
        end
        function varargout = m_solve(obj, varargin)
            %M_SOLVE - Accessing inv(M)
            %
            % See Also: HILUCSI4M_M_SOLVE, FACTORIZE
            [varargout{1:nargout}] = hilucsi4m_m_solve(obj.dbase, ...
                varargin{:});
        end
        function varargout = fgmres(obj, A, b, varargin)
            %FGMRES - Solving a system with FGMRES
            %
            % See Also: HILUCSI4M_FGMRES, FACTORIZE
            [varargout{1:nargout}] = hilucsi4m_fgmres(obj.dbase, A, b, ...
                varargin{:});
        end
        function varargout = gmres_null(obj, A, b, varargin)
            %GMRES_NULL - Solving for left null space
            %
            % See Also: HILUCSI4M_GMRES_NULL
            [varargout{1:nargout}] = hilucsi4m_gmres_null(obj.dbase, A, b, ...
                varargin{:});
        end
        function clear(obj)
            %CLEAR - Clear the internal storage
            hilucsi4m_clear(obj.dbase);
        end
        function delete(obj)
            % destructor
            hilucsi4m_finalize(obj.dbase);
        end
    end
end

