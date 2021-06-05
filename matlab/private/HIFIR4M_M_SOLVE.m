function flag = HIFIR4M_M_SOLVE
%HIFIR4M_M_SOLVE - Get the flag for accessing inv(M)

flag = int32(HIFIR4M_FACTORIZE + 1);
end