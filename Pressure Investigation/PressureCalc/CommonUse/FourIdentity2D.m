function [Jd, Js, Ja, Jm] = FourIdentity2D()
% FourIdentity2D.m
% 4th-order identity tensors
%
% Output: Jd, 2*2*2*2 matrix
%         Js, 2*2*2*2 matrix
%         Ja, 2*2*2*2 matrix
%         Jm, 2*2*2*2 matrix   
%--------------------------------------------------------------------------
%  delta
   delta = [1 0 ; 0 1];
  
%  allocate J, Jm
   J  = zeros(2,2,2,2);
   J1 = zeros(2,2,2,2);
   Jm = zeros(2,2,2,2);

%  Eqn(3) in Jiang(2014)
   for i = 1:2
       for j = 1:2
           for k = 1:2
               for l = 1:2
                   J(i,j,k,l)  = delta(i,k)*delta(j,l);
                   J1(i,j,k,l) = delta(j,k)*delta(i,l);
                   Jm(i,j,k,l) = delta(i,j)*delta(k,l)/2;
               end
           end
       end
   end
   Js = 0.5*(J + J1);
   Ja = 0.5*(J - J1);
   Jd = Js - Jm; 
   end