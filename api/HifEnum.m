classdef HifEnum
    % Enumeration data types of HIFIR4m
    properties (Constant)
        CREATE = int32(0);
        GET = int32(1);
        CLEAR = int32(2);
        CHECK = int32(3);
        DESTROY = int32(4);
        FACTORIZE = int32(5);
        M_SOLVE = int32(6);
        M_MULTIPLY = int32(7);
        EXPORT_DATA = int32(8);
        QUERY = int32(9);
    end
end
