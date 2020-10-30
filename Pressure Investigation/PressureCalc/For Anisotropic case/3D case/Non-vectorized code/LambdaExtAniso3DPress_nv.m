function p = LambdaExtAniso3DPress_nv(a,x,Cm,Alp,Bet,ww,psi,ww1)
 % This function is to calculate auxiliary tensor for pressure(Lambda) for an exterior point around
 % an inclusion embedded in an anisotropic viscous medium. 
 % Input: a is a 3*1 matrix, three semi-axes of the inclusion
 %        Cm is the stiffness tensor of the medium 
 %        x is the 3*1 matrix with external coordinates
 %        Alp and Bet are the Gaussian grid nodes for the product Gaussian
 %        quadrature (See Qu et al 2016)
 %        Alp is 1*N1^2 matrix containing Gaussian grid nodes corresponding to theta
 %        Bet is 1*N1^2 matrix containing Gaussian grid nodes corresponding to phi
 %        ww 1*N1^2 matrix containing Gaussian weights corresponding Alp and
 %        Bet,
 %        Here, N1 is the number of nodes and weights; The total grid nodes
 %        becomes N1^2.
 %        psi is N2*1 matrix containing Gaussian grid node corresponding to psi
 %        ww1 is N2*1 matrix containing Gaussian weights corresponding to psi
 %        Here, N2 is the number of nodes and weights;
 
 % Output: p is the Auxiliary tensor for Pressure Calculation (denoted by
 % 3*3 matrix).
 
   p = zeros(3,3);
    for i = 1:3
        for j = 1:3
            p(i,j) = AuxTensorPressOut_nv([i,j], a,x,0,2*pi, 0,pi,Alp, Bet, ww,Cm,psi,ww1);
        end
    end            
 
end