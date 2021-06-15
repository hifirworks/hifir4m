function is = crs_rowind(rowptr, js) %#codegen
% Obtain the row indices from the compressed row format.
%
%     is = crs_rowind(rowptr, js)
%
% See also crs_createFromAIJ

%#codegen -args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0),[inf,1])}

is = coder.nullcopy(zeros(size(js),class(js)));

nrows = int32(length(rowptr))-1;
for i=1:nrows
    for j = rowptr(i) : rowptr(i+1) - 1
        is(j) = i;
    end
end
