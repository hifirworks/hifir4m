
function bs = backsolve(R, bs, cend)
% Perform backward substitution.
%     bs = backsolve(R, bs)
%     bs = backsolve(R, bs, cend)
%     bs = backsolve(R, bs, cend, ws)
%  Compute bs = (triu(R(1:cend,1:cend))\bs) ./ ws;
%  The right-hand side vector bs can have multiple columns.
%  By default, cend is size(R,2), and ws is ones.

coder.inline('always');

if nargin<3
    cend = int32(size(R,2)); 
end

if isempty(coder.target)
    bs(1:cend,:) = R(1:cend, 1:cend) \ bs(1:cend,:);
else    
    for kk=1:int32(size(bs,2))
        for jj=cend:-1:1
            for ii=jj+1:cend
                bs(jj,kk) = bs(jj,kk) - R(jj,ii) * bs(ii,kk);
            end

            assert(R(jj,jj)~=0);
            bs(jj,kk) = bs(jj,kk) / R(jj,jj);
        end
    end
end
end
