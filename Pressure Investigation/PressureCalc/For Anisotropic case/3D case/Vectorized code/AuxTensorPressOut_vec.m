function T= AuxTensorPressOut_vec(sub, x,y, a, b, c, d,alp, bet, ww,Cm,psi,ww1)
                      
% Input: 
%        sub is a 1*2 matrix, 2 subscripts of the 2th order Auxiliary tensor for pressure (Lambda);
%        x is a 3*1 matrix, 3 semi-axes of an inclusion;
%        y is 3*1 matrix, 3 coordinates of the external point   
%        a,b is the integrating range of theta; c,d is the integrating range of phi;
%        Alp and Bet are the Gaussian grid nodes for the product Gaussian quadrature
%        Alp is 1*N1^2 matrix containing Gaussian grid nodes corresponding to theta
%        Bet is 1*N1^2 matrix containing Gaussian grid nodes corresponding to phi
%        ww 1*N1^2 matrix containing Gaussian weights corresponding Alp and Bet,
%        Here, N1 is the number of nodes and weights; The total grid nodes becomes (N1^2).
%        Cm is the stiffness tensor of the matrix
%        psi is N2*1 matrix containing Gaussian grid node corresponding to psi
%        ww1 is N2*1 matrix containing Gaussian weights corresponding to psi
%        Here, N2 is the number of nodes and weights;


% Output: T is a scalar value, corresponding to the elements of Auxiliary
% tensor(Lambda) for pressure.
 
 dp = DP_vec(x,y,alp,bet,Cm,psi,ww1);   % dp is a 3*3*(N1^2) matrix
 
 ft= squeeze(dp(sub(1),sub(2),:))'.*sin(bet);
     
 k   = ww.*ft;
 T = (d-c)*(b-a)*sum(k)/4;

end









