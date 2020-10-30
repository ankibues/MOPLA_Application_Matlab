function [S,p,PI] = SnpIn(x,Jd,Js,Ja)
 % This function is to calculate Eshelby tensors for Interior points
 % using integration
 
 % Input: x is a 3*1 matrix, three semi-axes of the inclusion
 
 % Output: S,p 
 % According to Jiang's mathcad sheet(EshelbyTensorIn&Out3,SPIP(x))
 
 % 31/1/2017 (AB) (used GGQ instead of APGQ)
  gp = 50;
 [Alp,Bet,ww] = GaussGGQ(gp);
   T = zeros(3,3,3,3);
   p = zeros(3,3);
    for i = 1:3
        for j = 1:3
            for m = 1:3
                for n = 1:3
                    T(i,j,m,n) = GGQ('fgin',[i,j,m,n],x,0,2*pi,0,pi,Alp,Bet,ww);
                end             
            end
            p(i,j) = GGQ('pgin',[i,j],x,0,2*pi,0,pi,Alp,Bet,ww);
        end
    end
                
   z  = Contract(T,2*Jd);
   S  = Contract(Js,z);
   PI = Contract(Ja,z);

end