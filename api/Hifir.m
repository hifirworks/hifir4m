classdef Hifir < handle
    % Handle for a HIFIR object

    properties
        hdl
        params
        A
    end

    methods
        function h = Hifir(hdl, params, A)
            % Constructor for Hifir object
            if nargin > 0
                h.hdl = hdl;
            else
                h.hdl= [];
            end

            if nargin > 1
                h.params= params;
            else
                h.params = '';
            end

            if nargin > 2
                h.A = A;
            else
                h.A = [];
            end
        end

        function varargout = apply(h, x, varargin)
            % Apply the preconditioner
            [varargout{1:nargout}] = hifApply(h, x, varargin{:});
        end

        function update(h, A)
            % Update coefficient matrix used in iterative refinement
            h = hifUpdate(h, A);
        end

        function varargout = refactorize(h, S, varargin)
            % Refactorization
            [varargout{1:nargout}] = hifRefactorize(h, S, varargin{:});
        end

        function delete(h)
            % Destructor for handle
            if ~isempty(h.hdl)
                hifir4m_mex(HifEnum.DESTROY, h.hdl);
            end
        end
    end
end
