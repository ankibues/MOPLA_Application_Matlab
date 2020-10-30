function T=  GreenTensorOutAniso_vec(sub, x,y, a, b, c, d,alp, bet, ww,Cm,psi,ww1)
                      
% for calculating the components of Green's interaction tensor T for
% exterior anisotropic case

% Input: 
%        sub is a 1*4 matrix, 4 subscripts of the 4th order tensor T;
%        x is a 3*1 matrix, 3 semi-axes of an inclusion;
%        y is 3*1 matrix, 3 coordinates of the external point   
%        a,b is the integrating range of theta; c,d is the integrating
%        range of phi;e,f is the integrating range of psi
%        Alp, Bet, and Gamma are Gaussian points for three dimensions
%        w is the weights
%        Cm is the stiffness tensor of the matrix

% Output: T is a scalar
 % theta= (pi)*alp + pi;   % changing the limits   1*n
 % phi= (pi/2)*bet + pi/2;
 
  dp = DT_vec(x,y,alp,bet,Cm,psi,ww1);   % 3*3*3*3*N1
 
  ft= squeeze(dp(sub(1),sub(2),sub(3),sub(4),:))'.*sin(bet);
     
  k   = ww.*ft;
  T = (d-c)*(b-a)*sum(k)/4;

end

%{
 function fr = funIn(sub,a,x,alp,bet,Cm)

end
%}










