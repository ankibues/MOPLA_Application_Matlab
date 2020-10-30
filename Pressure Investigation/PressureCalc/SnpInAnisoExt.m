function [S,p,PI] = SnpInAnisoExt(a,Jd,Ja,ang,Nn,m,x,Cm)
 % This function is to calculate Eshelby tensors for Interior points
 % using integration (Modify it later for calculating eshelby tensors)
 
 % Input: a is a 2*1 matrix, two semi-axes of the inclusion
 %        ang is angle phi, x is 2*n matrix where n is the number of external points 
 % Output: S,p 

 
 %  (AB) 
   S = zeros(2,2,2,2,numel(x)/2);
   PI= zeros(2,2,2,2,numel(x)/2);
   T = zeros(2,2,2,2,numel(x)/2);
   p = zeros(2,2,numel(x)/2);
   for o=1: (numel(x)/2)
       
       
       for i = 1:2
           for j = 1:2
               for m = 1:2
                   for n = 1:2
                       [T(i,j,m,n,o),~] = GGQAnisoExt([i,j,m,n],a,x(:,o),ang,Nn,m);
                   end             
               end
               [~,p(i,j,o)] = GGQAnisoExt([i,j,m,n],a,x(:,o), ang,Nn,m);
           end
       end
       z  = Contract2D(T(:,:,:,:,o),Cm);
       S(:,:,:,:,o)  = Contract2D(Jd,z);
       PI(:,:,:,:,o) = Contract2D(Ja,z);            
   end 

end