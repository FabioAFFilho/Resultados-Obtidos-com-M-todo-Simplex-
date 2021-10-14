function [x,y,iteracao,fo] = SimplexRevisado(N,b,c)
   [m,n] = size(N);
   #b = b + rand(m,1)*10e-12;
   b = sparse(b); #+ rand(m,1)*10e-12;
   
   bneg = find(b<0);
   N(bneg,:) = -N(bneg,:);
   b(bneg) = -b(bneg);
   s = find(b > -10e-8 && b < 10e-8);
   b(s) = 0;
   
   Au = [N eye(m)];
   Au = sparse(Au);
   d = [zeros(1,n) ones(1,m)];      
   d = sparse(d);
   nbase = [1:n];
   base = [n+1:n+m];
   [L,U,P,Q] = lu(Au(:,base));
              L = sparse(L);
              U = sparse(U);  
              P = sparse(P);
              Q = sparse(Q);
  xB = SolSistemLU(L, U, P*b);
  xB = Q*xB;
   #################Inicia o simplex revisado##################
  for k = 1:m+n
        

        y = SolSistemLU(U',L',Q'*d(base)');
        y = P'*y;
        s = find(y > -10e-8 && y < 10e-8);
        y(s) = 0;
        
        rN = d(nbase)' - Au(:,nbase)'*y;
        s = find(rN > -10e-8 && rN < 10e-8);
        rN(s) = 0;
        
        [val,j] = min(rN);

        if (val >= -10e-8);
           [x,y,iteracao,fo] = SimplexFaseII(N,b,c,base,nbase,xB,k);
          return
        endif

          Ntil = SolSistemLU(L,U,P*Au(:,nbase(j)));
          Ntil = Q*Ntil;
          
          s = find(Ntil > -10e-7 && Ntil < 10e-7);
          Ntil(s) = 0;
          
          if Ntil<= 10e-8;
             printf("Problema ilimitado!");
            return
          endif
          
          s=find(Ntil>10e-8);
          [val,i] = min(xB(s)./Ntil(s));
          

           aux10 = base(s(i));
           base(s(i)) = nbase(j);
           nbase(j)  = aux10;

           [L,U,P,Q] = lu(Au(:,base));
              L = sparse(L);
              U = sparse(U);  
              P = sparse(P);
              Q = sparse(Q);
           xB = SolSistemLU(L,U,P*b);
           xB = Q*xB;
           
           s=find(xB > -10e-8 && xB < 10e-8);
           xB(s) = 0;
           
           [valminxB,indice] = min(xB);
           if valminxB < -10e-8;
             printf("Valor infactível em xB")
             [val,i] = min(xB)
             return
           endif
       k=k+1;
     endfor
endfunction