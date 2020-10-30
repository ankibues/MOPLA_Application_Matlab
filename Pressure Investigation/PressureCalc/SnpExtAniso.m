function [S,p,PI] = SnpExtAniso(a,Jd,Ja,ang,Nn,m,x,Cm,Alp,ww)    
 % This function is to calculate Eshelby tensors for Interior points
 % using integration 
 
 % Input: a is a 2*1 matrix, two semi-axes of the inclusion
 % x is a 3*1 matrix with 3 coordinates
 
 % Output: S,p 

 
 %  (AB) 
  
 %[Alp,ww] = Gauss1(1000);
   
   
   
       T = zeros(2,2,2,2);
       p= zeros(2,2);
       X= [a,x];
       for i = 1:2
           for j = 1:2
               for k = 1:2
                   for l = 1:2
                       T(i,j,k,l) = GGQAnisoInandEx('fgex',[i,j,k,l],X,Alp,ww,ang,m,Nn);
                   end             
               end
           end
       end
       for ii= 1:2
           for jj=1:2
               p(ii,jj) = GGQAnisoInandEx('pgex',[ii,jj],X,Alp,ww,ang,m,Nn);
           end
       end
      
   z  = Contract2D(T,Cm);
   S  = Contract2D(Jd,z);
   PI = Contract2D(Ja,z);   
   
  
   end 
