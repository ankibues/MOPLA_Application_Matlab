function T = GGQS_PressCalc(sub,x,y, gp,Cm)
% Global Gaussian Quadrature on Sphere
% Numerical integration on the sphere according to Atkinson K.(1982)
% 2. Product Gaussian Quadrature
% range of integration: theta - 0,2*pi; phi - 0,pi;

% Input: f is a string, 'in' or 'ex';
%        sub is a 1*4 matrix, 4 subscripts of the 4th order tensor T;
%        x, if f is 'in', x is a 3*1 matrix, 3 semi-axes of an inclusion;
%           if f is 'ex', x is a 3*2 matrix, 1st column is 3 semi-axes of 
%                  an inclusion, 2nd column is an coordinates of an point;
%        gp is the number of gaussian points;

% Output: T is a scalar

  [p,w] = Gauss(gp);
  theta = (pi/gp):(pi/gp):(2*pi);
 
  phi = acos(p);  
  Alp = zeros(2*gp,gp);
  Bet = zeros(2*gp,gp);
  
  for i = 1:2*gp
    for j = 1:gp
        Alp(i,j) = theta(i);                                          
        Bet(i,j) = phi(j);
     end
  end
  
  w  = w';
  ww = repmat(w,1,2*gp);
      ft = funIn(sub,x,y,Alp,Bet,Cm);
  
  
  k  = ww.*ft;
  T  = (pi/gp)*sum(k);
end



function fr = funIn(sub,a,x,alpha,beta,Cm)

   theta = reshape(alpha,1,[]);
   phi   = reshape(beta,1,[]);
   dp = DPvectorized(a,x,theta,phi,Cm);   % 3*3*N1
   fr= squeeze(dp(sub(1),sub(2),:))'.*sin(phi);
   
end

