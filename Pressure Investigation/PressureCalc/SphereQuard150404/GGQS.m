function T = GGQS(f, sub, x, gp)
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
  
  if strcmp(f,'in')
      ft = funIn(sub,x,Alp,Bet);
  elseif strcmp(f,'ex')
      ft = funEx(sub,x(:,1),x(:,2),Alp,Bet);
  else
      fprintf('please choose between in and ex');
  end
  
  k  = ww.*ft;
  T  = (pi/gp)*sum(k);
end



function fr = funIn(sub,a,alpha,beta)

   theta = reshape(alpha,1,[]);
   phi   = reshape(beta,1,[]);
   n = numel(theta);
  
   xi = [cos(theta).*sin(phi); sin(theta).*sin(phi); cos(phi)];
   aaa = a(1)*a(2)*a(3)/(4*pi);
   
   aa = repmat(a,1,n);
   xia  = aa.*xi;
   xiasq = xia.^2;
   nxia = sqrt(sum(xiasq));
   rho = aaa./nxia.^3;
   
    if strcmp(ff,'fgin')
       del = Delta(sub(1),sub(3));
       delta = repmat(del,1,n);
       fr = rho.*xi(sub(4),:).*xi(sub(2),:).*(delta - xi(sub(1),:).*xi(sub(3),:));
   elseif strcmp(ff,'pgin')
       fr = -rho.*xi(sub(1),:).*xi(sub(2),:);
   else
       fprintf('need to choose a function (fgin or pgin)');
   end
   
   
end

function fr = funEx(sub,a,x,Alp,Bet)
  
  theta = reshape(Alp,1,[]);
  phi   = reshape(Bet,1,[]);
  n = numel(theta);
  
  xi = [cos(theta).*sin(phi); sin(theta).*sin(phi); cos(phi)];
  aaa     = a(1)*a(2)*a(3)/(4*pi);
  
  aa = repmat(a,1,n);
  xx = repmat(x,1,n);
  
  epsa  = aa.*xi;
  A     = xi./aa;
  r     = xx-epsa;
  
  rsq = r.^2;
  nr = sqrt(sum(rsq));
  
  temp = aaa./nr.^3;
  rho(1,:) = temp;
  rho(2,:) = temp;
  rho(3,:) = temp;
  rhoA = rho.*A;
  
  if strcmp(ff,'fgex')
      fr = 0.5*rhoA(sub(2),:).*(-Delta(sub(1),sub(3))*r(sub(4),:)...
               +Delta(sub(1),sub(4))*r(sub(3),:)+Delta(sub(3),sub(4))*r(sub(1),:)-...
               3*(r(sub(1),:).*r(sub(3),:).*r(sub(4),:))./(nr.^2));
  elseif strcmp(ff,'pgex')  
      fr = r(sub(1),:).*rhoA(sub(2),:);
  else
      fprintf('need to choose a function (fgex or pgex)');
  end
end

function re = Delta(m,n)
  % kroneckerDelta
  if m == n
      re = 1;
  else
      re = 0;
  end
end