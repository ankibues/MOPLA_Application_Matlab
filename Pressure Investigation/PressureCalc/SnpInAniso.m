function [S,p,PI] = SnpInAniso(x,Jd,Ja,ang,Nn,m,Cm,Alp,ww)
 % This function is to calculate Eshelby tensors for Interior points
 % using integration 
 
 % Input: x is a 2*1 matrix, two semi-axes of the inclusion
 
 % Output: S,p 

 
 %  (AB) 
  %[Alp,ww]=lgwt(1000,0,pi)
 
 %[Alp,ww] = Gauss1(8000);

   T = zeros(2,2,2,2);
   p = zeros(2,2);
    for i = 1:2
        for j = 1:2
            for k = 1:2
                for l = 1:2
                    T(i,j,k,l) = GGQAnisoInandEx('fgin',[i,j,k,l],x,Alp,ww,ang,m,Nn);
                end             
            end
        end
    end
   
    for ii=1:2
        for jj=1:2       
            p(ii,jj) = GGQAnisoInandEx('pgin',[ii,jj],x,Alp,ww,ang,m,Nn);
        end
    end
    z  = Contract2D(T,Cm);
    S  = Contract2D(Jd,z);
    PI = Contract2D(Ja,z);            
  

end