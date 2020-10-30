function T = T_ExtAniso3D_vec(a,x,Cm,Alp,Bet,ww,psi,ww1)
 % This function is to calculate Green's interaction tensor for exterior point in an
 % anisotropic case
 
 % using integration
 
 % Input: a is a 3*1 matrix, three semi-axes of the inclusion
 %  Cm stiffness tensor 
 % x is the 3*1 matrix with external coordinates
 
 % Output: T 


 
   T = zeros(3,3,3,3);
    for i = 1:3
        for j = 1:3
            for k=1:3
                for l=1:3
                    T(i,j,k,l) = GreenTensorOutAniso_vec([i,j,k,l], a,x,0,2*pi, 0,pi,Alp, Bet, ww,Cm,psi,ww1);
                end
            end
        end
    end            
 
end