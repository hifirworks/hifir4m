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
        KSP_SOLVE = int32(7);
        KSP_NULL = int32(8);
        EXPORT_DATA = int32(9);
        M_SOLVE2 = int32(10);
        M_MULTIPLY = int32(11);
        QUERY = int32(12);
    end
end
