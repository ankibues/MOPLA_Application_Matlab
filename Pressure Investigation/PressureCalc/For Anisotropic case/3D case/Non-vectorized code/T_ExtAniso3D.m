function T = T_ExtAniso3D(a,x,Cm,Alp,Bet,ww,psi,ww1)
 % This function is to calculate Green's interaction tensor for exterior point in an
 % anisotropic case
 % Note: this is for non-vectorized version, and really slow. Avoid using
 % it for general calculations. Only made for verification of vectorized
 % code.

 
 % Input: a is a 3*1 matrix, three semi-axes of the inclusion
 %  Cm stiffness tensor 
 % x is the 3*1 matrix with external coordinates
 
 % Output: T 


 
   T = zeros(3,3,3,3);
    for i = 1:3
        for j = 1:3
            for k=1:3
                for l=1:3
                    T(i,j,k,l) = GreenTensorOutAniso([i,j,k,l], a,x,0,2*pi, 0,pi,Alp, Bet, ww,Cm,psi,ww1);
                end
            end
        end
    end            
 
end