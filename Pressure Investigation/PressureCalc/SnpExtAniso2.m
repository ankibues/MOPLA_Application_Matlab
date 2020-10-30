function [S,p,PI] = SnpExtAniso2(a,Jd,Ja,ang,Nn,m,x,Cm)    
 % This function is to calculate Eshelby tensors for Interior points
 % using integration 
 
 % Input: a is a 2*1 matrix, two semi-axes of the inclusion
 
 % Output: S,p 

 
 %  (AB) 
  
tic
[Alp,ww] = Gauss1(2000);
toc   
S = zeros(2,2,2,2,numel(x)/2);
   PI= zeros(2,2,2,2,numel(x)/2);
   T = zeros(2,2,2,2,numel(x)/2);
   p = zeros(2,2,numel(x)/2);
   tic
   for o=1:(numel(x)/2)
       V= zeros(2,2,2,2);
       X1= x(:,o);
       X= [a,X1];
       for i = 1:2
           for j = 1:2
               for k = 1:2
                   for l = 1:2
                      V(i,j,k,l) = GGQAnisoInandExoptimized('fgex',[i,j,k,l],X,Alp,ww,ang,m,Nn);
                   end             
               end
           end
       end
       T(:,:,:,:,o)= V;
       
       VV= zeros(2,2);
       for ii= 1:2
           for jj=1:2
               VV(ii,jj) = GGQAnisoInandExoptimized('pgex',[ii,jj],X,Alp,ww,ang,m,Nn);
           end
       end
      p(:,:,o)= VV;
   z  = Contract2D(V,Cm);
   S(:,:,:,:,o)  = Contract2D(Jd,z);
   PI(:,:,:,:,o) = Contract2D(Ja,z);   
   end
   toc
   end 
