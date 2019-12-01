classdef HILUCSI
    %HILUCSI - OOP class representation of HILUCSI4M project
    %
    % Description:
    %   This class provides another way to use this package in an OOP
    %   fashion, which can be benefitial for certain programmers.
    
    properties (Access = protected)
        dbase
    end
    
    methods
        function obj = HILUCSI(varargin)
            %HILUCSI - Construct an instance with internal database
            %
            % Syntax:
            %   hl = HILUCSI
            %   hl = HILUCSI(true)
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
            flag = obj.dbase.is_mixed;
        end
        
        function varargout = factorize(obj, A, varargin)
            %FACTORIZE - Perform MLILU factorization
            %
            % See Also: HILUCSI4M_FACTORIZE, M_SOLVE
            [varargout{1:nargout}] = hilucsi4m_factorize(obj.dbase, A, ...
                varargin{:});
        end
        
        function varargout = m_solve(obj, b)
            %M_SOLVE - Accessing inv(M)
            %
            % See Also: HILUCSI4M_M_SOLVE, FACTORIZE
            [varargout{1:nargout}] = hilucsi4m_m_solve(obj.dbase, b);
        end
        
        function varargout = fgmres(obj, A, b, varargin)
            %FGMRES - Solving a system with FGMRES
            %
            % See Also: HILUCSI4M_FGMRES, FACTORIZE
            [varargout{1:nargout}] = hilucsi4m_fgmres(obj.dbase, A, b, ...
                varargin{:});
        end
        
        function delete(obj)
            % destructor
            hilucsi4m_finalize(obj.dbase);
        end
    end
end

