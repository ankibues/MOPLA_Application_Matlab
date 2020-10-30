function [S,PI,T] = SnpIn_all(x,Cm,Jd,Ja)
 % This function calculates Green tensor for all RDPs at once.
 
 % Input: x is a 3*N matrix, three semi-axes of N inclusions.
 
 % Output: S,p 
 % According to Jiang's mathcad sheet(EshelbyTensorIn&Out3,SPIP(x))
 
 % 31/1/2017 (AB) (used GGQ instead of APGQ)
  gp = 100;    % what should be the appropriate number of gp ?
  [Alp,Bet,ww] = GaussGGQ(gp);
  [~,N] = size(x);
   T = zeros(3,3,3,3,N);
   %p = zeros(3,3,N);
    for i = 1:3
        for j = 1:3
            for m = i:3
                for n = j:3
                    T(i,j,m,n,:) = GGQ_all2([i,j,m,n],x,0,2*pi,0,pi,Alp,Bet,ww,Cm);
                    T(m,j,i,n,:) = T(i,j,m,n,:);            
                    T(i,n,m,j,:) = T(i,j,m,n,:);
                    T(m,n,i,j,:) = T(i,j,m,n,:); 
                end             
            end
           % p(i,j) = GGQ('pgin',[i,j],x,0,2*pi,0,pi,Alp,Bet,ww);
        end
    end
                
   z  = Contract_vectorized(T,Cm);
   S  = Contract_vectorized2(Jd,z);
   PI = Contract_vectorized2(Ja,z);
   
end