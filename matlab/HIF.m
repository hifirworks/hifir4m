classdef HIF
    %HIF - OOP class representation of HIFIR4M project
    %
    % Description:
    %   This class provides another way to use this package in an OOP
    %   fashion, which can be benefitial for certain programmers.

    properties (Access = protected)
        dbase  % internal database structure
    end

    methods
        function obj = HIF(varargin)
            %HIF - Construct an instance with internal database
            %
            % Syntax:
            %   hl = HIF
            %   hl = HIF(true)
            %   hl = HIF(___, true)
            %
            % Description:
            %   hl = HIF simpily constructs an instance without mixed
            %   precision support.
            %
            %   hl = HIF(true) in constrast to the constructor above,
            %   this one enables mixed precision.
            %
            % See Also: HIFIR4M_INITIALIZE
            obj.dbase = hifir4m_initialize(varargin{:});
        end
        function flag = is_mixed(obj)
            %IS_MIXED - Check if the HIFIR instance is mixed precision
            %
            % See Also: ID
            flag = obj.dbase.is_mixed;
        end
        function flag = is_real(obj)
            %IS_REAL - Check if the HIFIR instance is real arithmetic
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
            flag = hifir4m_empty(obj.dbase);
        end
        function varargout = factorize(obj, A, varargin)
            %FACTORIZE - Perform MLILU factorization
            %
            % See Also: HIFIR4M_FACTORIZE, M_SOLVE
            [varargout{1:nargout}] = hifir4m_factorize(obj.dbase, A, ...
                varargin{:});
        end
        function hilu = export(obj, varargin)
            %EXPORT - Export internal data in C++ to MATLAB
            %
            % See Also: HIFIR4M_EXPORT
            hilu = hifir4m_export(obj.dbase, varargin{:});
        end
        function varargout = solve(obj, varargin)
            %SOLVE - Accessing inv(M)
            %
            % See Also: HIFIR4M_SOLVE, FACTORIZE
            [varargout{1:nargout}] = hifir4m_solve(obj.dbase, ...
                varargin{:});
        end
        function varargout = hifir(obj, varargin)
            %HIFIR - Accessing inv(M) with iterative refinement
            %
            % See Also: HIFIR4M_SOLVE, FACTORIZE
            [varargout{1:nargout}] = hifir4m_hifir(obj.dbase, ...
                varargin{:});
        end
        function varargout = gmres(obj, A, b, varargin)
            %GMRES - Solving a system with GMRES
            %
            % See Also: HIFIR4M_GMRES, FACTORIZE
            [varargout{1:nargout}] = hifir4m_gmres(obj.dbase, A, b, ...
                varargin{:});
        end
        function varargout = fgmres_null(obj, A, b, varargin)
            %FGMRES_NULL - Solving for null space
            %
            % See Also: HIFIR4M_FGMRES_NULL
            [varargout{1:nargout}] = hifir4m_fgmres_null(obj.dbase, A, b, ...
                varargin{:});
        end
        function clear(obj)
            %CLEAR - Clear the internal storage
            hifir4m_clear(obj.dbase);
        end
        function delete(obj)
            % destructor
            hifir4m_finalize(obj.dbase);
        end
    end
end

