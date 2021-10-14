function [x,y,iteracao,fo] = SimplexRevisado(N,b,c)
   N = sparse(N);
   [val,i] = min(b);
   [m,n] = size(N);
   b = b + rand(m,1)*10e-12;
   bneg = find(b<0);
   s = find(b > -10e-8 && b < 10e-8);
   b(s) = 0;
   N(bneg,:) = -N(bneg,:);
   b(bneg) = -b(bneg);
   B = eye(m);
   Au = [N B];
   d = [zeros(1,n) ones(1,m)];        
   nbase = [1:n];
   base = [n+1:n+m];
   [L,U,P,Q] = lu(Au(:,base));
   #xB = SolSistemLU(L,U,b,m);
   xB = Au(:,base)\b;
   Ntil = zeros(m,1);
   #################Inicia o simplex revisado##################
  for k = 1:m+n
        
        
        # Solução de y=Au(:,base)'\d(base) com fatoração (LU)' = U'L'
        y = SolSistemLU(U',L',Q'*d(base)',m);
        y = P'*y;
        #y=Au(:,base)'\d(base)';
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
        
          #entra na base        
          #Ntil = SolSistemLU(L,U,P*Au(:,nbase(j)),m);
          Ntil = SolSistemLU(L,U,P*Au(:,nbase(j)),m);
          Ntil = Q*Ntil;
          
          s = find(Ntil > -10e-7 && Ntil < 10e-7);
          Ntil(s) = 0;
          
          if max(Ntil)<= 10e-8;
             printf("Problema ilimitado!");
            return
          endif
          
          s=find(Ntil>0);
          [val,i] = min(xB(s)./Ntil(s));
          

           aux10 = base(s(i));
           base(s(i)) = nbase(j);
           nbase(j)  = aux10;

           [L,U,P,Q] = lu(Au(:,base));

           xB = SolSistemLU(L,U,P*b,m);
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