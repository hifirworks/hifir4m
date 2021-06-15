function s = vec_sqnorm2(a)
coder.inline('always');

if isempty(coder.target)
    s = a'*a;
else
    s = 0;
    for ii=1:int32(numel(a))
        s = s + conj(a(ii)) * a(ii);
    end
end
end
