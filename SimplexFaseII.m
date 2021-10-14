function [x,y,iteracao,fo] = SimplexFaseII(N,b,c,base,nbase,xB,k)
   [m,n] = size(N);
   N = sparse(N);
   nbase = [1:n];
   nbase(base) = [];
   x = zeros(n:1);
   Ntil = zeros(m,1);
   [L,U,P,Q] = lu(N(:,base));
   q = k
  for q = 1:n+m+k
        # Solução de y=N(:,base)'\c(base) com fatoração (LU)' = U'L'
        #y1 = SolSistemTransLU(U',L',P,c(base)',m);
        
        y = SolSistemLU(U',L',Q'*c(base)',m);
        y = P'*y;
        
        s = find(y > -10e-8 && y < 10e-8);
        y(s) = 0;
        
        rN = c(nbase)' - N(:,nbase)'*y;
        s = find(rN > -10e-8 && rN < 10e-8);
        rN(s) = 0;
        [val,j] = min(rN);
        
        if (val >= -10e-8)
          iteracao = q;
          fo = c(base)*xB;
          x(base) = xB;
          fprintf("Variável Basica fase II : [");
          fprintf("%g,", xB(1:end-1));
          fprintf("%g]\n", xB(end));
          fprintf("\n");
          fprintf("y =  [");
          fprintf("%g,", y(1:end-1));
          fprintf("%g]\n", y(end));
          fprintf("\n");
          fprintf("Numero iterações %d", q)
          fprintf("\n");
          fprintf("Valor funcao objetivo %d", c(base)*xB)
          fprintf("\n");
          fprintf("\n");
          return
        endif
        
          Ntil = SolSistemLU(L,U,P*N(:,nbase(j)),m);
          Ntil = Q*Ntil;
          
          s = find(Ntil > -10e-7 && Ntil < 10e-7);
          Ntil(s) = 0;
          s=find(Ntil>0);
          [val,i] = min(xB(s)./Ntil(s));
           if max(Ntil)<= 10e-8;
             printf("Problema ilimitado!");
            return
           endif
           
           aux10 = base(s(i));
           base(s(i)) = nbase(j);
           nbase(j)  = aux10;
           
           [L,U,P,Q] = lu(N(:,base));
           xB = SolSistemLU(L,U,P*b,m);
           xB = Q*xB;

            s = find(xB > -10e-8 && xB < 10e-8);
            xB(s) = 0;
            [valminxB,indice] = min(xB);
              if valminxB < -10e-8;
                printf("DEU PAU");
                fprintf("\n");
                fprintf("\n");
                fprintf("Base factível fase II : [");
                fprintf("%g,", xB(1:end-1));
                fprintf("%g]\n", xB(end));
                fprintf("\n");
                fprintf("valor de xB minimo %d", valminxB);
                fprintf("Numero iterações %d", q);
                fprintf("\n");
                fprintf("\n");
                return
             endif
           c(base)*xB;
       q=q+1;
     endfor
endfunction