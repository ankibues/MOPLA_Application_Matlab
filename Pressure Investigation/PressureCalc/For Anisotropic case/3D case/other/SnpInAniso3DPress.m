function p = SnpInAniso3DPress(x,Cm)
 % This function is to calculate Eshelby tensors for Interior points
 % For inclusion in an anisotropic medium
 
 % Input: x is a 3*1 matrix, three semi-axes of the inclusion
 
 % Output: S,p 
 % According to Jiang's mathcad sheet(EshelbyTensorIn&Out3,SPIP(x))
 
 % 31/1/2017 (AB) (used GGQ instead of APGQ)
  gp = 50;
 [Alp,Bet,ww] = GaussGGQ(gp);
   
   p = zeros(3,3);
    for i = 1:3
        for j = 1:3
            p(i,j) = GGQ3DPress('pgin',[i,j],x,0,2*pi,0,pi,Alp,Bet,ww,Cm);
        end
    end            
 
end