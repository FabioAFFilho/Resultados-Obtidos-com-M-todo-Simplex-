function [x,y,iteracao,fo] = SimplexFaseII(N,b,c,base,nbase,xB,k)
   [m,n] = size(N);
   N = sparse(N);
   b = sparse(b);
   c = sparse(c);
   x = zeros(1,n);

   [L,U,P,Q] = lu(N(:,base));
      L = sparse(L);
      U = sparse(U);  
      P = sparse(P);
      Q = sparse(Q);
   q = k;
   
  for q = 1:100*(n+m+k)
      nbase = [1:n];
      nbase(base) = [];
        #calculando y
        y = SolSistemLU(U',L',Q'*c(base)');
        y = P'*y;
        s = find(y > -10e-8 && y < 10e-8);
        y(s) = 0;

        #calculando o custo relativo rN
        rN = c(nbase)' - N(:,nbase)'*y;
        s = find(rN > -10e-8 && rN < 10e-8);
        rN(s) = 0;
        [val,j] = min(rN);
        x(base) = xB;
        #Verificando se é uma solução ótima
        if (val > -10e-8)
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
        
        #Encontrando Ntil
          Ntil = SolSistemLU(L,U,P*N(:,nbase(j)));
          Ntil = Q*Ntil;
          
          s = find(Ntil > -10e-7 && Ntil < 10e-7);
          Ntil(s) = 0;
          s=find(Ntil>10e-8);
          [val,i] = min(xB(s)./Ntil(s));
           if Ntil<= -10e-8;
             printf("Problema ilimitado!");
            return
           endif
           
           # Atualizando a Base
           aux10 = base(s(i));
           base(s(i)) = nbase(j);
           nbase(j)  = aux10;
           
           [L,U,P,Q] = lu(N(:,base));
              L = sparse(L);
              U = sparse(U);  
              P = sparse(P);
              Q = sparse(Q);
              
           # Atualizando xB
           xB = SolSistemLU(L,U,P*b);
           xB = Q*xB;
            
            s = find(xB > -10e-8 && xB < 10e-8);
            xB(s) = 0;
            [valminxB,indice] = min(xB);
              if valminxB < -10e-8;
                iteracao = q;
                fo = c(base)*xB;
                printf("Valor infactível");
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