function p = coder
% Function for backward compatibility with older version of MATLAB and Octave.

p.target = '';

p.extrinsic = @nullfunc;
p.cinclude = @nullfunc;
p.varsize = @nullfunc;
p.inline = @nullfunc;
p.allowpcode = @nullfunc;
p.cstructname = @nullfunc;
p.ismatlabthread = true;
p.unroll = @copy;
p.resize = @copy;
p.nullcopy = @copy;
p.ignoreConst = @copy;

p.opaque = @opaque;
p.typeof = @copy;

end

function nullfunc(varargin)
end

function A=copy(A, varargin)
end

function x=opaque(~, value)
if nargin == 2
  x = value;
else
  x = '';
end;
end
