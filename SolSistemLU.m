function [x] = SolSistemLU(L,U,b,m)             
        y=L\b;
       x = U\y;
       return
endfunction