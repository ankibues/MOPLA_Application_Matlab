function [S,p,PI] = SnpInAniso1(x,Jd,Js,Ja,ang,Nn,m,Cm)
 % This function is to calculate Eshelby tensors for Interior points
 % using integration 
 
 % Input: x is a 2*1 matrix, two semi-axes of the inclusion
 
 % Output: S,p 

 
 %  (AB)(this approach needs correction) 
 
   
   T = zeros(2,2,2,2);
   p = zeros(2,2);
    for i = 1:2
        for j = 1:2
            for m = 1:2
                for n = 1:2
                    [T(i,j,m,n),~] = GGQAniso1([i,j,m,n],x,ang,Nn,m);
                end             
            end
            [~,p(i,j)] = GGQAniso1([i,j,m,n],x,ang,Nn,m);
        end
    end
    z  = Contract2D(T,Cm);
    S  = Contract2D(Jd,z);
    PI = Contract2D(Ja,z);            
  

end